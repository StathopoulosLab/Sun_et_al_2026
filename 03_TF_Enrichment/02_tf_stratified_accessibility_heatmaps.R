if (file.exists("./paths_config.r")) source("./paths_config.r")  # optional project-local path overrides
################################################################################
# STEP 3 HEATMAPS — TF-STRATIFIED, PER-CONDITION DAR SETS
#
# PURPOSE:
#   Generates heatmaps for a given condition's DAR sets with an additional
#   stratification layer by TF binding status:
#     • Runt-bound  : DARs overlapping Runt ChIP-seq peaks
#     • Cic-bound   : DARs overlapping Cic ChIP-seq peaks (if available)
#     • Both-bound  : DARs overlapping both Runt and Cic
#     • Unbound     : DARs with no overlap to either TF
#
# STRUCTURE PER DAR SET:
#   {group}/{stem}/all_dars/           — all tracks, sequential colormap
#   {group}/{stem}/gained_open/        — all tracks, grpA genotype colors
#   {group}/{stem}/gained_close/       — all tracks, grpB genotype colors
#   {group}/{stem}/focal/              — only the 2 comparison tracks
#   {group}/{stem}/tf_stratified/
#     runt_bound/                      — Runt-bound DARs only
#     cic_bound/                       — Cic-bound DARs only
#     both_bound/                      — Runt+Cic co-bound DARs
#     unbound/                         — No TF overlap
#   {group}/{stem}/beds/               — BED files (including TF strata)
#
# SIGN CONVENTION:
#   NEGATE_LFC_STEMS lists the direct contrasts whose log2FC sign needs
#   flipping to match the numerator convention used in the figures.
#   All temporal/derived sets are correctly signed upstream — not negated.
#
# TF BED FILES:
#   Runt : ./data/chipseq/RuntAb_IgG_p0.05_reproducible.bed
#   Cic  : ./data/tf_partners/CicsfGFP_GSE130584_MACS2_peaks.bed
#
# CONFIG:
#   FLANK_BP sets the computeMatrix flanking window (+/- bp around the
#   reference point) used for every heatmap/profile in this run. The
#   CONDITION_* block below holds the genotype panel, colors, and DAR-set
#   stems for one condition (here: the ventralized BOTv/BOTCv split) — swap
#   that block out (or source it from an external file) to re-run the exact
#   same script structure against a different condition/genotype panel.
#   Set FLANK_BP <- 2000 for the +/-2kb heatmaps, or 5000 for +/-5kb; both
#   were used in the paper.
################################################################################

options(stringsAsFactors = FALSE)

FLANK_BP       <- 2000  # +/- bp flanking window for computeMatrix (2000 = 2kb; use 5000 for 5kb)
CONDITION_TAG  <- "Split2_Vent"                    # short tag used in output paths
CONDITION_NAME <- "Ventralized (BOTv / BOTCv)"     # human-readable label used in headers/titles

cat(sprintf("\n=== STEP 3 HEATMAPS -- TF-STRATIFIED: %s (+/-%dkb) ===\n",
            toupper(CONDITION_NAME), FLANK_BP / 1000))
cat("Timestamp:", format(Sys.time()), "\n\n")

# ── PATHS ──────────────────────────────────────────────────────────────────────
# dar_base_dir / temp_base_dir point at this condition's Step1 DAR outputs;
# point them at the equivalent non-ventralized output directories to re-run
# the same script for that condition instead.
dar_base_dir  <- "../Generate_fresh_counts/Output/SevenGeno_nc14b"
temp_base_dir <- "../Generate_fresh_counts/Output/Split2_Vent_temporal_v3"
figures_dir   <- sprintf("./%s/figures_%s_TF_Stratified_%dkb", CONDITION_TAG, CONDITION_TAG, FLANK_BP / 1000)
merged_bw_dir <- "./BigWig_Diagnostic_All/merged/methodB"

RUNT_BED <- "./data/chipseq/RuntAb_IgG_p0.05_reproducible.bed"
CIC_BED  <- "./data/tf_partners/CicsfGFP_GSE130584_MACS2_peaks.bed"
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
threads <- 8

# ── SMOOTHING ──────────────────────────────────────────────────────────────────
# NOTE: computeMatrix has NO matrix-level smoothing flag -- there is no
# --smoothLength for computeMatrix (that flag belongs to bamCoverage/
# bamCompare, used upstream at bigWig-generation time, which is a separate
# step from this one). The only lever computeMatrix itself gives you for
# reducing profile jaggedness is --binSize: coarser bins average more signal
# per point. For real smoothing control (loess/Savitzky-Golay), use
# plot_smoothed_profile.R on the resulting matrix.gz instead.
bin_size  <- 75         # computeMatrix --binSize (was 10; coarser = smoother)
plot_type <- "se"       # mean +/- SE band on profiles

# ── CONTRAST / SCALE SETTINGS ──────────────────────────────────────────────────
# ZMAX_PERCENTILE: upper quantile used to clamp each heatmap's colorscale
# (computed per-matrix at run time). Lower = more contrast (clips more tail).
# Heatmaps with few peaks get an even tighter clamp via ZMAX_PERCENTILE_SMALL,
# since a handful of strong/broad peaks can otherwise wash out the rest of the
# rows in a small set.
ZMAX_PERCENTILE       <- 90    # default colorscale clamp: 90th percentile
ZMAX_PERCENTILE_SMALL <- 80    # tighter clamp for small peak sets
ZMAX_SMALL_N          <- 100   # peak-count threshold for the tighter clamp
ZMIN_FIXED            <- 0     # ATAC accessibility signal always starts at 0

# ── PROJECT PALETTE ────────────────────────────────────────────────────────────
GENO_COLORS <- list(
  BOTv         = "#2ca02c",
  BOTv_late    = "#1a6b1a",
  BOTCv        = "#e377c2",
  BOTCv_late   = "#b5369a"
)

# TF stratum colors for profile tracks
TF_COLORS <- list(
  runt_bound  = "#9467bd",   # purple  (Runt family)
  cic_bound   = "#ff7f0e",   # orange  (Cic family)
  both_bound  = "#d62728",   # red     (convergence)
  unbound     = "#7f7f7f"    # grey
)

# ── BIGWIG CONFIGURATION ───────────────────────────────────────────────────────
bigwig_config <- list(
  list(file  = file.path(merged_bw_dir, "BOTv_nc14b_mean_bgSubCPM.bigWig"),
       label = "BOTv_nc14b",
       color = GENO_COLORS$BOTv,
       cmap_all = "Greens"),
  list(file  = file.path(merged_bw_dir, "BOTv_nc14late_mean_bgSubCPM.bigWig"),
       label = "BOTv_nc14late",
       color = GENO_COLORS$BOTv_late,
       cmap_all = "Greens"),
  list(file  = file.path(merged_bw_dir, "BOTCv_nc14b_mean_bgSubCPM.bigWig"),
       label = "BOTCv_nc14b",
       color = GENO_COLORS$BOTCv,
       cmap_all = "RdPu"),
  list(file  = file.path(merged_bw_dir, "BOTCv_nc14late_mean_bgSubCPM.bigWig"),
       label = "BOTCv_nc14late",
       color = GENO_COLORS$BOTCv_late,
       cmap_all = "RdPu")
)

# ── VALIDATE BIGWIGS ───────────────────────────────────────────────────────────
cat("--- BigWig validation ---\n")
active_bw <- Filter(function(bw) {
  ok <- file.exists(bw$file)
  cat(sprintf("  %s  %s\n", if (ok) "OK     " else "MISSING", bw$file))
  ok
}, bigwig_config)

cat(sprintf("\n  %d / %d BigWig files available\n\n",
            length(active_bw), length(bigwig_config)))
if (length(active_bw) == 0)
  stop("No BigWig files found at: ", merged_bw_dir)

# ── VALIDATE TF BED FILES ──────────────────────────────────────────────────────
cat("--- TF BED validation ---\n")
use_runt <- file.exists(RUNT_BED)
use_cic  <- file.exists(CIC_BED)
cat(sprintf("  %s  Runt BED : %s\n", if (use_runt) "OK     " else "MISSING", RUNT_BED))
cat(sprintf("  %s  Cic BED  : %s\n", if (use_cic)  "OK     " else "MISSING", CIC_BED))
if (!use_runt && !use_cic) {
  stop("Neither Runt nor Cic BED files found. Cannot stratify.\n",
       "Expected:\n  ", RUNT_BED, "\n  ", CIC_BED)
}
cat("\n")

# ── FULL-PANEL STRINGS ─────────────────────────────────────────────────────────
bw_files_str     <- paste(sapply(active_bw, `[[`, "file"),    collapse = " ")
bw_labels_str    <- paste(sapply(active_bw, `[[`, "label"),   collapse = " ")
bw_colors_str    <- paste(paste0("'", sapply(active_bw, `[[`, "color"), "'"), collapse = " ")
bw_cmaps_all_str <- paste(sapply(active_bw, `[[`, "cmap_all"), collapse = " ")
bw_colorlist_str <- paste(
  sapply(active_bw, function(bw) paste0("'white,", bw$color, "'")),
  collapse = " ")

cat("Active tracks:\n")
for (bw in active_bw)
  cat(sprintf("  %-22s  %s  (cmap_all=%s)\n", bw$label, bw$color, bw$cmap_all))
cat("\n")

# ── DAR SET DEFINITIONS ────────────────────────────────────────────────────────
dar_sets <- list(

  # Direct comparisons (nc14b)
  list(stem = "BOTv_vs_BOTCv",
       dir = dar_base_dir, group = "direct",
       focal_bw_labels = c("BOTv_nc14b", "BOTCv_nc14b")),

  # Direct comparisons (nc14late)
  list(stem = "BOTv_vs_BOTCv_nc14late",
       dir = temp_base_dir, group = "direct",
       focal_bw_labels = c("BOTv_nc14late", "BOTCv_nc14late")),

  # Background / tolrm9 sets
  list(stem = "tolrm9_inferred_bg",
       dir = dar_base_dir, group = "background",
       focal_bw_labels = c("BOTv_nc14b", "BOTCv_nc14b")),
  list(stem = "tolrm9_empirical_bg",
       dir = dar_base_dir, group = "background",
       focal_bw_labels = c("BOTv_nc14b", "BOTCv_nc14b")),
  list(stem = "tolrm9_validated_both_methods",
       dir = dar_base_dir, group = "background",
       focal_bw_labels = c("BOTv_nc14b", "BOTCv_nc14b")),

  # Specificity sets
  list(stem = "BOTCv_above_inferred_bg",
       dir = dar_base_dir, group = "specificity",
       focal_bw_labels = c("BOTv_nc14b", "BOTCv_nc14b")),
  list(stem = "BOTCv_above_empirical_bg",
       dir = dar_base_dir, group = "specificity",
       focal_bw_labels = c("BOTv_nc14b", "BOTCv_nc14b")),
  list(stem = "BOTv_above_inferred_bg",
       dir = dar_base_dir, group = "specificity",
       focal_bw_labels = c("BOTv_nc14b", "BOTCv_nc14b")),
  list(stem = "BOTCv_only_crossvalidated",
       dir = dar_base_dir, group = "specificity",
       focal_bw_labels = c("BOTv_nc14b", "BOTCv_nc14b")),
  list(stem = "BOTv_only_crossvalidated",
       dir = dar_base_dir, group = "specificity",
       focal_bw_labels = c("BOTv_nc14b", "BOTCv_nc14b")),
  list(stem = "BOTv_BOTCv_divergent_vs_BOTR",
       dir = dar_base_dir, group = "specificity",
       focal_bw_labels = c("BOTv_nc14b", "BOTCv_nc14b")),

  # Temporal sets
  list(stem = "BOTv_temporal",
       dir = temp_base_dir, group = "temporal",
       focal_bw_labels = c("BOTv_nc14b", "BOTv_nc14late")),
  list(stem = "BOTCv_temporal",
       dir = temp_base_dir, group = "temporal",
       focal_bw_labels = c("BOTCv_nc14b", "BOTCv_nc14late")),
  list(stem = "BOTv_BOTCv_stable_both_timepoints",
       dir = temp_base_dir, group = "temporal",
       focal_bw_labels = NULL),
  list(stem = "BOTv_BOTCv_emerging_late_only",
       dir = temp_base_dir, group = "temporal",
       focal_bw_labels = c("BOTv_nc14late", "BOTCv_nc14late")),
  list(stem = "BOTv_BOTCv_resolving_early_only",
       dir = temp_base_dir, group = "temporal",
       focal_bw_labels = c("BOTv_nc14b", "BOTCv_nc14b")),
  list(stem = "temporal_divergent_BOTv_vs_BOTCv",
       dir = temp_base_dir, group = "temporal",
       focal_bw_labels = NULL)
)

# Only these two direct contrasts need sign flip (BOTv-numerator convention)
NEGATE_LFC_STEMS <- c(
  "BOTv_vs_BOTCv",
  "BOTv_vs_BOTCv_nc14late"
)

# ── I/O RETRY HELPERS ───────────────────────────────────────────────────────────
# Wrap filesystem-touching operations (wc -l, file.copy) with a short retry
# loop. On a syncing/networked working directory (e.g. a cloud-synced folder),
# file creation/reads can transiently time out mid-batch. A raw system() call
# that times out returns character(0) instead of a count, which used to crash
# as.integer()/if(n > 0) outright. These wrappers retry a few times with a
# short pause, and return NA (rather than erroring) if genuinely stuck, so one
# transient hiccup doesn't kill an hours-long batch run. Defined up-front so
# BED validation (which runs immediately below) can use them too.
RETRY_ATTEMPTS  <- 4
RETRY_PAUSE_SEC <- 2

safe_wc_l <- function(path, attempts = RETRY_ATTEMPTS, pause = RETRY_PAUSE_SEC) {
  for (i in seq_len(attempts)) {
    if (file.exists(path)) {
      out <- suppressWarnings(system(paste("wc -l <", shQuote(path)), intern = TRUE))
      if (length(out) > 0) {
        n <- suppressWarnings(as.integer(out[1]))
        if (!is.na(n)) return(n)
      }
    }
    if (i < attempts) {
      cat(sprintf("    [retry %d/%d] wc -l timed out/empty on %s, retrying...\n",
                  i, attempts, basename(path)))
      Sys.sleep(pause)
    }
  }
  cat(sprintf("    [WARN] could not count lines in %s after %d attempts -- treating as empty\n",
              basename(path), attempts))
  NA_integer_
}

safe_file_copy <- function(from, to, attempts = RETRY_ATTEMPTS, pause = RETRY_PAUSE_SEC) {
  for (i in seq_len(attempts)) {
    ok <- suppressWarnings(tryCatch(file.copy(from, to, overwrite = TRUE),
                                     error = function(e) FALSE))
    if (isTRUE(ok) && file.exists(to)) return(TRUE)
    if (i < attempts) {
      cat(sprintf("    [retry %d/%d] file.copy timed out writing %s, retrying...\n",
                  i, attempts, basename(to)))
      Sys.sleep(pause)
    }
  }
  cat(sprintf("    [WARN] failed to write %s after %d attempts -- skipping\n",
              basename(to), attempts))
  FALSE
}

# ── BED VALIDATION ─────────────────────────────────────────────────────────────
cat("--- DAR BED validation ---\n")
for (ds in dar_sets) {
  bed <- file.path(ds$dir, paste0(ds$stem, "_DARs.bed"))
  if (file.exists(bed)) {
    n <- safe_wc_l(bed)
    cat(sprintf("  OK      %-48s  %d regions\n", basename(bed), n))
  } else {
    cat(sprintf("  MISSING %-48s  (will skip)\n", basename(bed)))
  }
}
cat("\n")

# ── BED SPLITTING HELPER ───────────────────────────────────────────────────────
# Identical logic to Split2 v2. Priority: col4 > col5 > annotated.txt direction > log2FC.
split_dar_bed <- function(bed_path, stem, negate_lfc = FALSE) {
  df <- tryCatch(
    read.table(bed_path, header = FALSE, sep = "\t",
               stringsAsFactors = FALSE, fill = TRUE, comment.char = ""),
    error = function(e) { cat("  ERROR reading BED:", e$message, "\n"); NULL }
  )
  if (is.null(df) || nrow(df) == 0)
    return(list(open = NULL, close = NULL, method = "none"))

  if (negate_lfc)
    cat(sprintf("  %s: applying sign flip (BOTv-numerator convention)\n", stem))

  # Method 1: col4 direction string
  if (ncol(df) >= 4) {
    col4          <- tolower(as.character(df[, 4]))
    raw_open_idx  <- grepl("gained.open|gain_open|opening|_open|gain(?!ed.close)",
                           col4, perl = TRUE)
    raw_close_idx <- grepl("gained.close|gain_close|closing|_close|loss", col4)
    if (sum(raw_open_idx) + sum(raw_close_idx) > 0) {
      open_idx  <- if (negate_lfc) raw_close_idx else raw_open_idx
      close_idx <- if (negate_lfc) raw_open_idx  else raw_close_idx
      cat(sprintf("  %s: col4 direction%s -- open=%d  close=%d\n",
                  stem, if (negate_lfc) " [swapped]" else "",
                  sum(open_idx), sum(close_idx)))
      return(list(open = if (sum(open_idx) > 0) df[open_idx, ] else NULL,
                  close = if (sum(close_idx) > 0) df[close_idx, ] else NULL,
                  method = "col4"))
    }
  }

  # Method 2: col5 signed score
  if (ncol(df) >= 5) {
    scores <- suppressWarnings(as.numeric(df[, 5]))
    if (!all(is.na(scores)) &&
        any(scores > 0, na.rm = TRUE) && any(scores < 0, na.rm = TRUE)) {
      display_scores <- if (negate_lfc) -scores else scores
      open_idx  <- !is.na(display_scores) & display_scores > 0
      close_idx <- !is.na(display_scores) & display_scores < 0
      cat(sprintf("  %s: col5 %s score -- open=%d  close=%d\n",
                  stem, if (negate_lfc) "negated" else "signed",
                  sum(open_idx), sum(close_idx)))
      return(list(open = if (sum(open_idx) > 0) df[open_idx, ] else NULL,
                  close = if (sum(close_idx) > 0) df[close_idx, ] else NULL,
                  method = "col5"))
    }
  }

  # Method 3: _DARs_annotated.txt fallback
  ann_path <- sub("_DARs\\.bed$", "_DARs_annotated.txt", bed_path)
  if (file.exists(ann_path)) {
    ann <- tryCatch(
      read.table(ann_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                 quote = "\"", fill = TRUE, comment.char = "", row.names = NULL),
      error = function(e) NULL)

    # Reorder ann to match df's peak order BEFORE extracting direction/lfc,
    # rather than trusting that matching row counts implies matching row
    # order. Tries name match, then a normalized (_DAR token stripped) name
    # match -- some contrasts' annotated files insert this token and others
    # don't -- then coordinate match. Falls back to the original row-count-
    # only assumption only if none of those are possible (no usable name/
    # coordinate columns), same safety bar as before, just no longer the
    # ONLY way this can succeed.
    if (!is.null(ann) && nrow(ann) > 0) {
      bed_name <- if (ncol(df) >= 4) as.character(df[, 4]) else NULL
      ann_name_col <- intersect(c("name","peak_name","peakID","peak_id","peak","ID"),
                                colnames(ann))[1]
      reordered <- FALSE
      if (!is.na(ann_name_col) && !is.null(bed_name)) {
        idx <- match(bed_name, as.character(ann[[ann_name_col]]))
        if (sum(!is.na(idx)) >= 0.5*nrow(df)) {
          ann <- ann[idx, , drop = FALSE]; reordered <- TRUE
          cat(sprintf("  %s: annotated file reordered via name match (%d/%d)\n",
                     stem, sum(!is.na(idx)), nrow(df)))
        } else {
          norm <- function(x) gsub("_DAR(?=_|$)", "", x, perl=TRUE)
          idx_n <- match(norm(bed_name), norm(as.character(ann[[ann_name_col]])))
          if (sum(!is.na(idx_n)) >= 0.5*nrow(df)) {
            ann <- ann[idx_n, , drop = FALSE]; reordered <- TRUE
            cat(sprintf("  %s: annotated file reordered via normalized name match, _DAR token stripped (%d/%d)\n",
                       stem, sum(!is.na(idx_n)), nrow(df)))
          }
        }
      }
      if (!reordered) {
        chr_col <- intersect(c("chr","chrom","seqnames"), colnames(ann))[1]
        st_col  <- intersect(c("start"), colnames(ann))[1]
        en_col  <- intersect(c("end"), colnames(ann))[1]
        if (!is.na(chr_col) && !is.na(st_col) && !is.na(en_col)) {
          bed_key <- paste(df[,1], df[,2], df[,3], sep=":")
          ann_key <- paste(ann[[chr_col]], ann[[st_col]], ann[[en_col]], sep=":")
          idx_c <- match(bed_key, ann_key)
          if (sum(!is.na(idx_c)) >= 0.5*nrow(df)) {
            ann <- ann[idx_c, , drop = FALSE]; reordered <- TRUE
            cat(sprintf("  %s: annotated file reordered via genomic coordinates (%d/%d)\n",
                       stem, sum(!is.na(idx_c)), nrow(df)))
          }
        }
      }
      # Original fallback: only trust position if row counts genuinely agree
      # and none of the above reordering was possible.
      ann_usable <- reordered || nrow(ann) == nrow(df)
    } else {
      ann_usable <- FALSE
    }

    if (ann_usable) {
      dir_col <- intersect(c("direction", "Direction"), colnames(ann))[1]
      if (!is.na(dir_col)) {
        dirs <- tolower(trimws(ann[[dir_col]]))
        raw_open_idx  <- !is.na(dirs) & dirs == "gained-open"
        raw_close_idx <- !is.na(dirs) & dirs == "gained-close"
        if (sum(raw_open_idx) + sum(raw_close_idx) > 0) {
          open_idx  <- if (negate_lfc) raw_close_idx else raw_open_idx
          close_idx <- if (negate_lfc) raw_open_idx  else raw_close_idx
          cat(sprintf("  %s: annotated direction%s -- open=%d  close=%d\n",
                      stem, if (negate_lfc) " [swapped]" else "",
                      sum(open_idx), sum(close_idx)))
          return(list(open = if (sum(open_idx) > 0) df[open_idx, ] else NULL,
                      close = if (sum(close_idx) > 0) df[close_idx, ] else NULL,
                      method = "annotated_direction"))
        }
      }

      lfc_col <- intersect(c("log2FC", "log2FoldChange", "logFC", "LFC"),
                           colnames(ann))[1]
      if (!is.na(lfc_col)) {
        lfc         <- suppressWarnings(as.numeric(ann[[lfc_col]]))
        display_lfc <- if (negate_lfc) -lfc else lfc
        if (!all(is.na(display_lfc)) &&
            any(display_lfc > 0, na.rm = TRUE) &&
            any(display_lfc < 0, na.rm = TRUE)) {
          open_idx  <- !is.na(display_lfc) & display_lfc > 0
          close_idx <- !is.na(display_lfc) & display_lfc < 0
          cat(sprintf("  %s: annotated log2FC%s -- open=%d  close=%d\n",
                      stem, if (negate_lfc) " [negated]" else "",
                      sum(open_idx), sum(close_idx)))
          return(list(open = if (sum(open_idx) > 0) df[open_idx, ] else NULL,
                      close = if (sum(close_idx) > 0) df[close_idx, ] else NULL,
                      method = "annotated_lfc"))
        }
      }
    }
  }

  cat(sprintf("  %s: no direction info -- only all_dars / TF strata generated\n", stem))
  list(open = NULL, close = NULL, method = "none")
}

# ── TF STRATIFICATION HELPER ───────────────────────────────────────────────────
# Uses bedtools intersect to split a DAR BED into up to four TF-binding strata.
# Writes BED files and returns named list of paths (NULL if stratum is empty).
# Requires: bedtools on PATH.
stratify_by_tf <- function(bed_path, stem, bed_dir) {
  if (is.null(bed_path) || !file.exists(bed_path)) return(NULL)

  n_total <- safe_wc_l(bed_path)
  if (is.na(n_total) || n_total == 0) return(NULL)

  strata <- list()

  safe_intersect <- function(a, b, invert = FALSE, out_path = NULL) {
    flag <- if (invert) " -v" else ""
    # -u reports each feature in A once if any overlap with B
    cmd  <- sprintf("bedtools intersect -a %s -b %s%s -u > %s",
                    shQuote(a), shQuote(b), flag, shQuote(out_path))
    system(cmd)
    n <- safe_wc_l(out_path)
    if (!is.na(n) && n > 0) out_path else NULL
  }

  # Runt-bound (possibly also Cic-bound — we keep the union for the runt_ track)
  if (use_runt) {
    p <- file.path(bed_dir, paste0(stem, "_runt_bound.bed"))
    strata$runt_bound <- safe_intersect(bed_path, RUNT_BED, out_path = p)
    if (!is.null(strata$runt_bound)) {
      n <- safe_wc_l(p)
      cat(sprintf("    [TF] runt_bound: %d regions\n", n))
    }
  }

  # Cic-bound
  if (use_cic) {
    p <- file.path(bed_dir, paste0(stem, "_cic_bound.bed"))
    strata$cic_bound <- safe_intersect(bed_path, CIC_BED, out_path = p)
    if (!is.null(strata$cic_bound)) {
      n <- safe_wc_l(p)
      cat(sprintf("    [TF] cic_bound: %d regions\n", n))
    }
  }

  # Both-bound (Runt AND Cic)
  if (use_runt && use_cic) {
    tmp  <- tempfile(fileext = ".bed")
    p    <- file.path(bed_dir, paste0(stem, "_both_bound.bed"))
    # First intersect with Runt, then with Cic
    cmd1 <- sprintf("bedtools intersect -a %s -b %s -u > %s",
                    shQuote(bed_path), shQuote(RUNT_BED), shQuote(tmp))
    system(cmd1)
    cmd2 <- sprintf("bedtools intersect -a %s -b %s -u > %s",
                    shQuote(tmp), shQuote(CIC_BED), shQuote(p))
    system(cmd2)
    unlink(tmp)
    n <- safe_wc_l(p)
    if (!is.na(n) && n > 0) {
      strata$both_bound <- p
      cat(sprintf("    [TF] both_bound: %d regions\n", n))
    }
  }

  # Unbound (no overlap with Runt or Cic)
  {
    p <- file.path(bed_dir, paste0(stem, "_unbound.bed"))
    tf_beds_used <- c(
      if (use_runt) RUNT_BED else NULL,
      if (use_cic)  CIC_BED  else NULL
    )
    # Chain: subtract each TF BED in turn
    tmp_in <- bed_path
    for (tf_bed in tf_beds_used) {
      tmp_out <- tempfile(fileext = ".bed")
      cmd <- sprintf("bedtools intersect -a %s -b %s -v > %s",
                     shQuote(tmp_in), shQuote(tf_bed), shQuote(tmp_out))
      system(cmd)
      if (tmp_in != bed_path) unlink(tmp_in)
      tmp_in <- tmp_out
    }
    copy_ok <- safe_file_copy(tmp_in, p)
    unlink(tmp_in)
    n <- if (copy_ok) safe_wc_l(p) else NA_integer_
    if (!is.na(n) && n > 0) {
      strata$unbound <- p
      cat(sprintf("    [TF] unbound:    %d regions\n", n))
    }
  }

  strata
}

# Returns the zMax percentile to use for a given BED — tighter clamp for
# small peak sets (n < ZMAX_SMALL_N), since a few strong/broad peaks can
# otherwise dominate the colorscale and wash out the rest of the rows.
zmax_pct <- function(bed_path) {
  if (is.null(bed_path) || !file.exists(bed_path)) return(ZMAX_PERCENTILE)
  n <- safe_wc_l(bed_path)
  if (is.na(n) || n < ZMAX_SMALL_N) ZMAX_PERCENTILE_SMALL else ZMAX_PERCENTILE
}

# ── COMMAND BLOCK HELPERS ──────────────────────────────────────────────────────

# Full-panel, colormap-per-track (all_dars)
make_block_cmap <- function(bed_path, tag, cmaps_str, odir,
                             files_str, labels_str, colors_str) {
  if (is.null(bed_path) || !file.exists(bed_path)) return(NULL)
  mat  <- file.path(odir, paste0(tag, "_matrix.gz"))
  heat <- file.path(odir, paste0(tag, "_heatmap.pdf"))
  prof <- file.path(odir, paste0(tag, "_profile.pdf"))
  c(
    paste0("# -- ", tag, " --"),
    paste("computeMatrix reference-point",
          sprintf("--referencePoint center -b %d -a %d", FLANK_BP, FLANK_BP),
          paste("--binSize", bin_size, "--skipZeros"),
          "-S", files_str,
          "-R", bed_path,
          "-o", mat,
          "--numberOfProcessors", threads),
    paste0("ZMAX_CURRENT=$(zmax_from_matrix ", mat, " ", zmax_pct(bed_path), ")"),
    paste("plotHeatmap",
          "-m", mat, "-o", heat,
          "--samplesLabel", labels_str,
          "--colorMap", cmaps_str,
          "--whatToShow 'heatmap and colorbar'",
          "--plotTitle", shQuote(tag),
          "--sortRegions descend",
          paste0("--zMin ", ZMIN_FIXED, " --zMax $ZMAX_CURRENT"),
          "--heatmapHeight 15 --heatmapWidth 4 --dpi 300"),
    paste("plotProfile",
          "-m", mat, "-o", prof,
          "--samplesLabel", labels_str,
          "--colors", colors_str,
          "--plotTitle", shQuote(tag),
          "--plotType", plot_type,
          "--perGroup --plotHeight 7 --plotWidth 8 --legendLocation upper-right"),
    ""
  )
}

# Directional block — genotype hex colors via --colorList
make_block_colorlist <- function(bed_path, tag, colorlist_str, odir,
                                  files_str, labels_str, colors_str) {
  if (is.null(bed_path) || !file.exists(bed_path)) return(NULL)
  mat  <- file.path(odir, paste0(tag, "_matrix.gz"))
  heat <- file.path(odir, paste0(tag, "_heatmap.pdf"))
  prof <- file.path(odir, paste0(tag, "_profile.pdf"))
  c(
    paste0("# -- ", tag, " --"),
    paste("computeMatrix reference-point",
          sprintf("--referencePoint center -b %d -a %d", FLANK_BP, FLANK_BP),
          paste("--binSize", bin_size, "--skipZeros"),
          "-S", files_str,
          "-R", bed_path,
          "-o", mat,
          "--numberOfProcessors", threads),
    paste0("ZMAX_CURRENT=$(zmax_from_matrix ", mat, " ", zmax_pct(bed_path), ")"),
    paste("plotHeatmap",
          "-m", mat, "-o", heat,
          "--samplesLabel", labels_str,
          "--colorList", colorlist_str,
          "--whatToShow 'heatmap and colorbar'",
          "--plotTitle", shQuote(tag),
          "--sortRegions descend",
          paste0("--zMin ", ZMIN_FIXED, " --zMax $ZMAX_CURRENT"),
          "--heatmapHeight 15 --heatmapWidth 4 --dpi 300"),
    paste("plotProfile",
          "-m", mat, "-o", prof,
          "--samplesLabel", labels_str,
          "--colors", colors_str,
          "--plotTitle", shQuote(tag),
          "--plotType", plot_type,
          "--perGroup --plotHeight 7 --plotWidth 8 --legendLocation upper-right"),
    ""
  )
}

# TF-stratum block: shared colormap (YlOrRd) across strata for easy comparison
make_block_tf_stratum <- function(bed_path, tag, odir,
                                   files_str, labels_str, colors_str,
                                   colorlist_str) {
  if (is.null(bed_path) || !file.exists(bed_path)) return(NULL)
  n <- safe_wc_l(bed_path)
  if (is.na(n) || n == 0) return(NULL)
  mat  <- file.path(odir, paste0(tag, "_matrix.gz"))
  heat <- file.path(odir, paste0(tag, "_heatmap.pdf"))
  prof <- file.path(odir, paste0(tag, "_profile.pdf"))
  c(
    paste0("# -- ", tag, " --"),
    paste("computeMatrix reference-point",
          sprintf("--referencePoint center -b %d -a %d", FLANK_BP, FLANK_BP),
          paste("--binSize", bin_size, "--skipZeros"),
          "-S", files_str,
          "-R", bed_path,
          "-o", mat,
          "--numberOfProcessors", threads),
    paste0("ZMAX_CURRENT=$(zmax_from_matrix ", mat, " ",
           if (n < ZMAX_SMALL_N) ZMAX_PERCENTILE_SMALL else ZMAX_PERCENTILE, ")"),
    # colorList keeps the genotype identity visible within each TF stratum
    paste("plotHeatmap",
          "-m", mat, "-o", heat,
          "--samplesLabel", labels_str,
          "--colorList", colorlist_str,
          "--whatToShow 'heatmap and colorbar'",
          "--plotTitle", shQuote(tag),
          "--sortRegions descend",
          paste0("--zMin ", ZMIN_FIXED, " --zMax $ZMAX_CURRENT"),
          "--heatmapHeight 15 --heatmapWidth 4 --dpi 300"),
    paste("plotProfile",
          "-m", mat, "-o", prof,
          "--samplesLabel", labels_str,
          "--colors", colors_str,
          "--plotTitle", shQuote(tag),
          "--plotType", plot_type,
          "--perGroup --plotHeight 7 --plotWidth 8 --legendLocation upper-right"),
    ""
  )
}

# ── GENERATE SHELL SCRIPT ──────────────────────────────────────────────────────
all_cmds <- c(
  "#!/usr/bin/env bash",
  "set -euo pipefail",
  "",
  "# ── per-matrix zMax computation ─────────────────────────────────────────────",
  "# zmax_from_matrix MATRIX.gz PERCENTILE",
  "# Reads the deepTools matrix, extracts all numeric signal values, and returns",
  "# the requested percentile. Used to clamp each heatmap's colorscale.",
  "zmax_from_matrix() {",
  "  local mat=$1 pct=${2:-90}",
  "  python3 - <<PYEOF",
  "import gzip, numpy as np",
  "vals = []",
  "with gzip.open('$mat', 'rt') as f:",
  "    for line in f:",
  "        if line.startswith('@'): continue",
  "        vals += [float(x) for x in line.split()[6:] if x != 'nan']",
  "print(round(np.percentile(vals, $pct), 2))",
  "PYEOF",
  "}",
  "",
  sprintf("# deepTools -- %s: TF-stratified heatmaps", CONDITION_NAME),
  paste0("# Generated: ", format(Sys.time())),
  paste0("# Tracks: ", bw_labels_str),
  "# Stratification: Runt-bound / Cic-bound / Both-bound / Unbound",
  "# Color convention:",
  "#   all_dars     -- genotype sequential colormaps (Greens / RdPu)",
  "#   gained_open  -- genotype hex colors via --colorList",
  "#   gained_close -- genotype hex colors via --colorList",
  "#   focal/       -- only the two comparison-relevant tracks",
  "#   tf_stratified/{runt,cic,both,unbound}/ -- genotype colorList within TF stratum",
  "#   Colorscale: zMax clamped per-heatmap at the 90th percentile (80th for",
  "#   peak sets < 100 regions) of matrix signal; zMin fixed at 0.",
  ""
)

n_written <- 0
n_split   <- 0

for (ds in dar_sets) {
  bed_all <- file.path(ds$dir, paste0(ds$stem, "_DARs.bed"))
  if (!file.exists(bed_all)) next

  cat(sprintf("\n--- Processing: %s  [%s] ---\n", ds$stem, ds$group))

  odir_root    <- file.path(figures_dir, ds$group, ds$stem)
  odir_all     <- file.path(odir_root, "all_dars")
  odir_open    <- file.path(odir_root, "gained_open")
  odir_close   <- file.path(odir_root, "gained_close")
  odir_focal   <- file.path(odir_root, "focal")
  odir_tf      <- file.path(odir_root, "tf_stratified")
  odir_tf_runt <- file.path(odir_tf, "runt_bound")
  odir_tf_cic  <- file.path(odir_tf, "cic_bound")
  odir_tf_both <- file.path(odir_tf, "both_bound")
  odir_tf_unb  <- file.path(odir_tf, "unbound")
  bed_dir      <- file.path(odir_root, "beds")

  for (d in c(odir_all, odir_open, odir_close, odir_focal,
              odir_tf_runt, odir_tf_cic, odir_tf_both, odir_tf_unb, bed_dir))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

  # -- Split by direction
  split         <- split_dar_bed(bed_all, ds$stem,
                                 negate_lfc = ds$stem %in% NEGATE_LFC_STEMS)
  bed_open_path  <- NULL
  bed_close_path <- NULL

  if (!is.null(split$open) && nrow(split$open) > 0) {
    bed_open_path <- file.path(bed_dir, "gained_open.bed")
    write.table(split$open, bed_open_path, sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    cat(sprintf("    saved gained_open.bed  (%d regions)\n", nrow(split$open)))
  }
  if (!is.null(split$close) && nrow(split$close) > 0) {
    bed_close_path <- file.path(bed_dir, "gained_close.bed")
    write.table(split$close, bed_close_path, sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    cat(sprintf("    saved gained_close.bed (%d regions)\n", nrow(split$close)))
  }
  if (split$method != "none") n_split <- n_split + 1

  # -- TF stratification on the full DAR set (all directions)
  cat(sprintf("  TF stratification of %s...\n", basename(bed_all)))
  strata <- stratify_by_tf(bed_all, ds$stem, bed_dir)

  # -- Focal track strings
  focal_labels <- ds$focal_bw_labels
  focal_bw <- if (!is.null(focal_labels)) {
    Filter(function(bw) bw$label %in% focal_labels, active_bw)
  } else {
    active_bw
  }
  focal_files_str     <- paste(sapply(focal_bw, `[[`, "file"),    collapse = " ")
  focal_labels_str    <- paste(sapply(focal_bw, `[[`, "label"),   collapse = " ")
  focal_colors_str    <- paste(paste0("'", sapply(focal_bw, `[[`, "color"), "'"), collapse = " ")
  focal_cmap_all_str  <- paste(sapply(focal_bw, `[[`, "cmap_all"), collapse = " ")
  focal_colorlist_str <- paste(
    sapply(focal_bw, function(bw) paste0("'white,", bw$color, "'")),
    collapse = " ")

  # -- Build shell commands
  all_cmds <- c(all_cmds,
    paste0("# ", strrep("=", 60)),
    paste0("# DAR SET: ", ds$stem, "  [", ds$group, "]"),
    paste0("# ", strrep("=", 60)),
    "",

    # 1. All tracks — sequential colormap (unchanged from Split2 v2)
    make_block_cmap(bed_all,
                    paste0(ds$stem, "_all_DARs"),
                    bw_cmaps_all_str, odir_all,
                    bw_files_str, bw_labels_str, bw_colors_str),

    # 2. Gained open — all tracks, genotype colors
    make_block_colorlist(bed_open_path,
                         paste0(ds$stem, "_gained_open"),
                         bw_colorlist_str, odir_open,
                         bw_files_str, bw_labels_str, bw_colors_str),

    # 3. Gained close — all tracks, genotype colors
    make_block_colorlist(bed_close_path,
                         paste0(ds$stem, "_gained_close"),
                         bw_colorlist_str, odir_close,
                         bw_files_str, bw_labels_str, bw_colors_str),

    # 4. Focal — all DARs, comparison tracks only
    make_block_cmap(bed_all,
                    paste0(ds$stem, "_focal_all"),
                    focal_cmap_all_str, odir_focal,
                    focal_files_str, focal_labels_str, focal_colors_str),

    # 5. Focal gained open
    make_block_colorlist(bed_open_path,
                         paste0(ds$stem, "_focal_open"),
                         focal_colorlist_str, odir_focal,
                         focal_files_str, focal_labels_str, focal_colors_str),

    # 6. Focal gained close
    make_block_colorlist(bed_close_path,
                         paste0(ds$stem, "_focal_close"),
                         focal_colorlist_str, odir_focal,
                         focal_files_str, focal_labels_str, focal_colors_str),

    # ── TF-STRATIFIED BLOCKS (new) ────────────────────────────────────────────
    "# -- TF stratified --",

    # 7. Runt-bound
    make_block_tf_stratum(
      strata$runt_bound,
      paste0(ds$stem, "_runt_bound"),
      odir_tf_runt, bw_files_str, bw_labels_str, bw_colors_str, bw_colorlist_str),

    # 8. Cic-bound
    make_block_tf_stratum(
      strata$cic_bound,
      paste0(ds$stem, "_cic_bound"),
      odir_tf_cic, bw_files_str, bw_labels_str, bw_colors_str, bw_colorlist_str),

    # 9. Both-bound (Runt + Cic)
    make_block_tf_stratum(
      strata$both_bound,
      paste0(ds$stem, "_both_bound"),
      odir_tf_both, bw_files_str, bw_labels_str, bw_colors_str, bw_colorlist_str),

    # 10. Unbound
    make_block_tf_stratum(
      strata$unbound,
      paste0(ds$stem, "_unbound"),
      odir_tf_unb, bw_files_str, bw_labels_str, bw_colors_str, bw_colorlist_str)
  )

  n_written <- n_written + 1
}

all_cmds    <- c(all_cmds, sprintf("echo 'All TF-stratified %s heatmaps done.'", CONDITION_TAG))
script_path <- file.path(figures_dir, sprintf("run_deeptools_%s_TF_Stratified_%dkb.sh", CONDITION_TAG, FLANK_BP / 1000))
writeLines(all_cmds, script_path)
Sys.chmod(script_path, "0755")

cat("\n")
cat("======================================================================\n")
cat(sprintf("STEP 3 COMPLETE -- TF-STRATIFIED: %s\n", toupper(CONDITION_NAME)))
cat("======================================================================\n\n")
cat(sprintf("Shell script:         %s\n\n", script_path))
cat(sprintf("DAR sets included   : %d / %d\n", n_written, length(dar_sets)))
cat(sprintf("  with open/close   : %d\n", n_split))
cat(sprintf("  all_dars only     : %d\n\n", n_written - n_split))
cat(sprintf("BigWig tracks       : %d / %d available\n\n",
            length(active_bw), length(bigwig_config)))
cat("Output per DAR set:\n")
cat("  {group}/{stem}/all_dars/                  -- all 4 tracks, sequential colormap\n")
cat("  {group}/{stem}/gained_open/               -- all 4 tracks, genotype hex colors\n")
cat("  {group}/{stem}/gained_close/              -- all 4 tracks, genotype hex colors\n")
cat("  {group}/{stem}/focal/                     -- only the 2 comparison tracks\n")
cat("  {group}/{stem}/tf_stratified/runt_bound/  -- Runt ChIP-overlapping DARs\n")
cat("  {group}/{stem}/tf_stratified/cic_bound/   -- Cic ChIP-overlapping DARs\n")
cat("  {group}/{stem}/tf_stratified/both_bound/  -- Runt+Cic co-bound DARs\n")
cat("  {group}/{stem}/tf_stratified/unbound/     -- No TF overlap DARs\n")
cat("  {group}/{stem}/beds/                      -- all BED files\n\n")
cat("Sign convention:\n")
cat("  NEGATE_LFC_STEMS = BOTv_vs_BOTCv, BOTv_vs_BOTCv_nc14late only.\n")
cat("  Temporal/derived sets NOT negated.\n\n")
cat("To run deepTools:\n")
cat(sprintf("  bash %s\n\n", script_path))
