################################################################################
# 02_peak_fate_classification_strict.r  —  MAXIMALLY STRICT PEAK-FATE CLASSIFICATION
#
# PURPOSE
# -------
# 01_peak_fate_classification.r's vent (BOTv vs BOTCv) classification does not rely
# solely on "is this locus a significant DAR at both nc14b and nc14late" --
# because nc14late's genotype-comparison test is far less powered than
# nc14b's (1,371 vs 6,323 significant calls, despite the two timepoints
# being ~30 min apart), v16 recovers additional loci through several
# explicit, separately-labeled fallback tiers (A3b/A3c/A5 in v16):
#
#   1. Unfiltered-trend fallback (A3b): an nc14b DAR with no significant
#      nc14late DAR match is still counted as "matched" if the RAW
#      (non-significance-filtered) nc14late log2FC clears an absolute
#      magnitude floor (|log2FC|>0.5), OR sits within 70-130% of the
#      nc14b log2FC in the same direction ("stability-consistency" rescue).
#   2. Trajectory-imputation fallback (A3c): loci still unmatched after (1)
#      get an nc14late genotype-difference ESTIMATED by combining each
#      genotype's own independently significant temporal DAR calls
#      (predicted = lfc_nc14b + BOTv's own shift - BOTCv's own shift),
#      again admitted either by an absolute floor (>0.6) or the same
#      stability-consistency band.
#   3. Emerged (de novo) support tiers (A5): beyond a reproducible
#      MACS2-called peak, v16 also accepts a single-replicate "union"
#      call, or a depth+magnitude floor (baseMean>=10, |log2FC|>=0.2), as
#      evidence a claimed Emerged peak is real and not near-zero-count
#      DESeq2 sign noise.
#   4. One HCR-validated locus added back below the standard ATAC
#      significance/magnitude bar entirely on independent orthogonal
#      evidence (hcr_validated_loci in v16).
#
# None of this is a bug -- each tier is real, principled, explicitly
# labeled per-locus (match_confidence / emerged_direction_support /
# validation_source), and documented in v16's own comments. But it does
# mean the class totals reported from v16's output are LARGER than what a
# reviewer would get by re-deriving "differentially accessible" from the
# significance-filtered DAR calls alone -- which is worth showing directly.
#
# This script does two things:
#
#   PART 1 (STRICT RECLASSIFICATION): re-derives the vent fate table using
#   ONLY loci with a directly significant DAR call at BOTH timepoints
#   (tier-1 join, same reciprocal-overlap logic as v16's A3) for
#   Maintained/Reversed/Converged, and ONLY reproducible-called-peak-
#   supported de novo calls for Emerged -- i.e. every fallback tier above
#   is switched off. No unfiltered-trend fallback, no trajectory
#   imputation, no union/depth rescue, no HCR allowlist addback.
#
#   PART 2 (ATTRIBUTION): if v16's own peak_fate_data.csv is available, it
#   is read back in and every locus is bucketed by WHICH tier it actually
#   needed (straight from its recorded match_confidence /
#   emerged_direction_support / validation_source), so the default-vs-
#   strict gap can be decomposed rather than just reported as a single
#   number.
#
# OUTPUTS (./Overview_Plots/peak_fate/BOTv_vs_BOTCv/STRICT/):
#   peak_fate_data.csv                  -- strict-only locus table. Named to
#                                           match v16's export exactly (not
#                                           suffixed) so 04_peak_fate_alluvial_plots_strict.r
#                                           can read it with zero changes --
#                                           living in STRICT/ is what
#                                           disambiguates it.
#   fate_comparison_default_vs_strict.csv
#   fate_comparison_default_vs_strict.txt
#   fate_comparison_barplot.pdf         -- default vs strict counts/fate
#
# Run 04_peak_fate_alluvial_plots_strict.r afterwards for a strict-mode alluvial (same
# plots as 03_peak_fate_alluvial_plots.r, reading/writing entirely within this
# STRICT/ folder instead of the default one).
#
# This script only touches PART A (vent). v16's Part B (non-vent temporal
# fate) has no equivalent fallback tiers, so it isn't affected by this
# distinction and doesn't need a strict counterpart.
#
# USAGE: run from the same working directory as 01_peak_fate_classification.r (same
# ./ATAC_DARS input layout). Run 01_peak_fate_classification.r first if you want
# PART 2's attribution breakdown -- PART 1 works standalone either way.
################################################################################

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)
})

################################################################################
# CONFIG -- keep in sync with 01_peak_fate_classification.r's CONFIG block
################################################################################

dar_dir  <- "./ATAC_DARS"
out_dir  <- "./Overview_Plots"

stem_nc14b         <- "BOTv_vs_BOTCv"
stem_nc14late      <- "BOTv_vs_BOTCv_nc14late"
grpA <- "BOTv";  lbl_grpA <- "BOT.v-biased"
grpB <- "BOTCv"; lbl_grpB <- "BOTC.v-biased"

# fate_split labels -- must match 01_peak_fate_classification.r / plot_alluvial_*.r
# exactly, so 04_peak_fate_alluvial_plots_strict.r can read this script's output directly.
lbl_mnt_A <- "Maintained BOT.v bias"     # Maintained, ends grpA-open
lbl_mnt_B <- "Maintained BOTC.v bias"    # Maintained, ends grpB-open
lbl_rev_A <- "Reversed to BOTC.v"        # Reversed, started grpA-open
lbl_rev_B <- "Reversed to BOT.v"         # Reversed, started grpB-open

min_recip_frac <- 0.25   # identical to v16 -- this is DAR-to-DAR overlap
                         # geometry, not a significance rescue, so it's kept

# Called-peak support, STRICT tier only: reproducible (multi-rep) calls.
# The union (single-rep) tier and the depth+magnitude rescue are
# deliberately NOT used here -- see header.
called_peak_dir <- "./Merged_Peak_BEDs"
called_peak_file <- function(geno, timepoint) {
  if (is.null(called_peak_dir)) return(NA_character_)
  file.path(called_peak_dir, sprintf("%s_%s_merged_reproducible.bed", geno, timepoint))
}

# Explicit construct/driver-locus exclusion -- kept identical to v16. This
# is an artifact exclusion, not a significance rescue, so strict mode still
# applies it (it can only ever REMOVE loci, never add them back).
emerged_exclude_loci <- GRanges(
  seqnames = c("chr2R"),
  ranges   = IRanges(start = c(17756000), end = c(17758500))
)

out_contrast        <- file.path(out_dir, "peak_fate", "BOTv_vs_BOTCv")
out_strict      <- file.path(out_contrast, "STRICT")
dir.create(out_strict, recursive = TRUE, showWarnings = FALSE)

################################################################################
# HELPERS (subset of v16's, unchanged -- copy here so this script is
# runnable standalone without sourcing v16)
################################################################################

read_dar_full <- function(stem, search_dirs = dar_dir) {
  for (d in search_dirs) {
    bed_path <- file.path(d, paste0(stem, "_DARs.bed"))
    if (!file.exists(bed_path)) next
    bed <- tryCatch(read.table(bed_path, sep = "\t", header = FALSE,
                                stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(bed) || nrow(bed) == 0) next
    n_bed <- nrow(bed)
    lfc   <- rep(NA_real_, n_bed)
    padj  <- rep(NA_real_, n_bed)
    depth <- rep(NA_real_, n_bed)

    for (ext in c("_DARs_annotated.txt", "_DARs_annotated_peaks.csv")) {
      ann_path <- file.path(d, paste0(stem, ext))
      if (!file.exists(ann_path)) next
      sep_ch <- if (grepl("\t", readLines(ann_path, 1, warn = FALSE), fixed = TRUE)) "\t" else ","
      ann <- tryCatch(read.table(ann_path, header = TRUE, sep = sep_ch,
                                  stringsAsFactors = FALSE, quote = "\"",
                                  fill = TRUE, comment.char = "", row.names = NULL),
                       error = function(e) NULL)
      if (is.null(ann)) next
      colnames(ann)[1] <- gsub("^[^[:alnum:]]+", "", colnames(ann)[1])
      lc <- intersect(c("log2FoldChange", "log2FC", "LFC"), colnames(ann))[1]
      pc <- intersect(c("padj", "adj.P.Val", "FDR", "p_adj"), colnames(ann))[1]
      dc <- intersect(c("baseMeaAveExpr", "baseMean", "AveExpr", "baseMean.AveExpr"), colnames(ann))[1]
      if (is.na(lc)) next

      bed_name <- if (ncol(bed) >= 4) as.character(bed[[4]]) else NULL
      id_col <- intersect(c("name", "peak_name", "peakID", "peak_id", "peak", "ID", "Row.names"),
                           colnames(ann))[1]

      matched <- FALSE
      if (!is.na(id_col) && !is.null(bed_name)) {
        idx <- match(bed_name, as.character(ann[[id_col]]))
        hit <- !is.na(idx)
        if (sum(hit) >= 0.5 * n_bed) {
          lfc[hit] <- as.numeric(ann[[lc]])[idx[hit]]
          if (!is.na(pc)) padj[hit] <- as.numeric(ann[[pc]])[idx[hit]]
          if (!is.na(dc)) depth[hit] <- as.numeric(ann[[dc]])[idx[hit]]
          matched <- TRUE
        } else {
          norm <- function(x) gsub("_DAR(?=_|$)", "", x, perl = TRUE)
          idx_n <- match(norm(bed_name), norm(as.character(ann[[id_col]])))
          hit_n <- !is.na(idx_n)
          if (sum(hit_n) >= 0.5 * n_bed) {
            lfc[hit_n] <- as.numeric(ann[[lc]])[idx_n[hit_n]]
            if (!is.na(pc)) padj[hit_n] <- as.numeric(ann[[pc]])[idx_n[hit_n]]
            if (!is.na(dc)) depth[hit_n] <- as.numeric(ann[[dc]])[idx_n[hit_n]]
            matched <- TRUE
          }
        }
      }
      if (!matched) {
        chr_col   <- intersect(c("chr", "chrom", "seqnames", "Chr", "Chrom"), colnames(ann))[1]
        start_col <- intersect(c("start", "Start", "chromStart"), colnames(ann))[1]
        end_col   <- intersect(c("end", "End", "chromEnd"), colnames(ann))[1]
        if (!is.na(chr_col) && !is.na(start_col) && !is.na(end_col)) {
          bed_key <- paste(bed[[1]], bed[[2]], bed[[3]], sep = ":")
          ann_key <- paste(ann[[chr_col]], ann[[start_col]], ann[[end_col]], sep = ":")
          idx <- match(bed_key, ann_key)
          hit <- !is.na(idx)
          if (sum(hit) >= 0.5 * n_bed) {
            lfc[hit] <- as.numeric(ann[[lc]])[idx[hit]]
            if (!is.na(pc)) padj[hit] <- as.numeric(ann[[pc]])[idx[hit]]
            if (!is.na(dc)) depth[hit] <- as.numeric(ann[[dc]])[idx[hit]]
            matched <- TRUE
          }
        }
      }
      if (!matched && nrow(ann) == n_bed) {
        lfc <- as.numeric(ann[[lc]])
        if (!is.na(pc)) padj <- as.numeric(ann[[pc]])
        if (!is.na(dc)) depth <- as.numeric(ann[[dc]])
        matched <- TRUE
      }
      if (matched) break
    }

    if (all(is.na(lfc)) && ncol(bed) >= 5) {
      raw5 <- suppressWarnings(as.numeric(bed[, 5]))
      if (any(raw5 < 0, na.rm = TRUE)) lfc <- raw5
    }
    out <- data.frame(
      chr = as.character(bed[[1]]), start = as.integer(bed[[2]]), end = as.integer(bed[[3]]),
      name = if (ncol(bed) >= 4) as.character(bed[[4]]) else paste0(stem, "_", seq_len(nrow(bed))),
      log2FC = as.numeric(lfc), padj = as.numeric(padj), depth = as.numeric(depth),
      stringsAsFactors = FALSE)
    return(out)
  }
  message("  [missing] ", stem); NULL
}

bed_to_gr <- function(df) {
  GRanges(seqnames = df$chr, ranges = IRanges(start = df$start + 1L, end = df$end),
          name = df$name, log2FC = df$log2FC, padj = df$padj,
          depth = if ("depth" %in% names(df)) df$depth else NA_real_)
}

best_recip_join <- function(gr_query, gr_subject, min_recip = 0.25) {
  hits <- findOverlaps(gr_query, gr_subject, minoverlap = 1L)
  if (length(hits) == 0) return(data.frame(i_q = integer(), i_s = integer(), recip_frac = numeric()))
  ov_w  <- width(pintersect(gr_query[queryHits(hits)], gr_subject[subjectHits(hits)]))
  recip <- pmin(ov_w / width(gr_query[queryHits(hits)]), ov_w / width(gr_subject[subjectHits(hits)]))
  hdf <- data.frame(i_q = queryHits(hits), i_s = subjectHits(hits), recip_frac = recip)
  hdf <- hdf[hdf$recip_frac >= min_recip, ]
  if (nrow(hdf) == 0) return(hdf)
  hdf %>% group_by(i_q) %>% slice_max(recip_frac, n = 1, with_ties = FALSE) %>%
    ungroup() %>% as.data.frame()
}

load_called_peaks <- function(geno, timepoint) {
  fp <- called_peak_file(geno, timepoint)
  if (is.na(fp) || !file.exists(fp)) return(NULL)
  pk <- tryCatch(read.table(fp, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                             comment.char = "#"), error = function(e) NULL)
  if (is.null(pk) || nrow(pk) == 0) return(NULL)
  tryCatch(GRanges(seqnames = pk[[1]], ranges = IRanges(pk[[2]] + 1L, pk[[3]])),
           error = function(e) NULL)
}

################################################################################
# PART 1 -- STRICT RECLASSIFICATION
################################################################################

cat(strrep("=", 70), "\nPART 1: STRICT-ONLY vent peak fate (no trend/imputation/union/depth/HCR rescue)\n", strrep("=", 70), "\n\n", sep = "")

df_nc14b    <- read_dar_full(stem_nc14b)
df_nc14late <- read_dar_full(stem_nc14late)
if (is.null(df_nc14b) || is.null(df_nc14late))
  stop("Cannot find required BED files for vent comparison")

cat(sprintf("  nc14b genotype DARs   : %d\n", nrow(df_nc14b)))
cat(sprintf("  nc14late genotype DARs: %d\n", nrow(df_nc14late)))

df_nc14b$direction_nc14b <- ifelse(df_nc14b$log2FC > 0, lbl_grpA, lbl_grpB)
df_nc14b$direction_nc14b[is.na(df_nc14b$log2FC)] <- "Unknown"

gr_nc14b    <- bed_to_gr(df_nc14b)
gr_nc14late <- bed_to_gr(df_nc14late)

# ── Tier-1-only join: significant DAR <-> significant DAR reciprocal overlap.
# This is IDENTICAL to v16's A3 first pass, before A3b/A3c run. Nothing else
# is applied.
join_geno <- best_recip_join(gr_nc14b, gr_nc14late, min_recip = min_recip_frac)
cat(sprintf("  Matched nc14b->nc14late (significant DAR both timepoints): %d / %d (%.1f%%)\n",
            nrow(join_geno), nrow(df_nc14b), 100 * nrow(join_geno) / nrow(df_nc14b)))

fate_df <- as.data.frame(df_nc14b[, c("chr", "start", "end", "name", "log2FC", "padj", "direction_nc14b")],
                          stringsAsFactors = FALSE)
colnames(fate_df)[colnames(fate_df) == "log2FC"] <- "lfc_nc14b"
colnames(fate_df)[colnames(fate_df) == "padj"]   <- "padj_nc14b"
fate_df$lfc_nc14late  <- NA_real_
fate_df$padj_nc14late <- NA_real_
fate_df$matched_nc14late <- FALSE
if (nrow(join_geno) > 0) {
  fate_df$lfc_nc14late[join_geno$i_q]     <- mcols(gr_nc14late[join_geno$i_s])$log2FC
  fate_df$padj_nc14late[join_geno$i_q]    <- mcols(gr_nc14late[join_geno$i_s])$padj
  fate_df$matched_nc14late[join_geno$i_q] <- TRUE
}

fate_df <- fate_df %>% mutate(
  fate = case_when(
    !matched_nc14late                                       ~ "Converged",
    matched_nc14late & is.na(lfc_nc14late)                  ~ "Converged",
    matched_nc14late & !is.na(lfc_nc14late) & sign(lfc_nc14late) == sign(lfc_nc14b) ~ "Maintained",
    matched_nc14late & !is.na(lfc_nc14late) & sign(lfc_nc14late) != sign(lfc_nc14b) ~ "Reversed",
    TRUE ~ "Converged"),
  direction_nc14late = case_when(
    fate == "Converged"                     ~ "Converged",
    !is.na(lfc_nc14late) & lfc_nc14late > 0  ~ lbl_grpA,
    !is.na(lfc_nc14late) & lfc_nc14late < 0  ~ lbl_grpB,
    TRUE                                     ~ "Converged"))

cat("\n  Strict fate breakdown (Maintained/Reversed/Converged only, pre-Emerged):\n")
print(table(fate_df$direction_nc14b, fate_df$fate))

# ── Emerged, strict tier only: reproducible called-peak support required.
# No union (single-rep) tier, no depth+magnitude rescue, no HCR allowlist.
emerged_mask <- !overlapsAny(gr_nc14late, gr_nc14b, minoverlap = 1L)
n_emerged <- sum(emerged_mask)
cat(sprintf("\n  De novo peaks at nc14late (pre-support-check): %d\n", n_emerged))

if (n_emerged > 0) {
  gr_em <- gr_nc14late[emerged_mask]
  claimed_direction <- ifelse(!is.na(mcols(gr_em)$log2FC) & mcols(gr_em)$log2FC > 0, grpA, grpB)

  peaks_grpA <- load_called_peaks(grpA, "nc14late")
  peaks_grpB <- load_called_peaks(grpB, "nc14late")
  has_called_peaks <- !is.null(peaks_grpA) || !is.null(peaks_grpB)

  if (!has_called_peaks) {
    cat("  [NOTE] No reproducible called-peak files found -- strict mode cannot verify\n",
        "  Emerged calls against called peaks and will report them WITHOUT the peak-\n",
        "  support filter applied (i.e. this run is not actually as strict as it could\n",
        "  be for the Emerged category specifically). Provide the same\n",
        "  ./Merged_Peak_BEDs/<geno>_nc14late_merged_reproducible.bed files v16 uses to\n",
        "  enable this check.\n")
    supported <- rep(TRUE, n_emerged)
  } else {
    supported <- vapply(seq_len(n_emerged), function(i) {
      pk <- if (claimed_direction[i] == grpA) peaks_grpA else peaks_grpB
      !is.null(pk) && overlapsAny(gr_em[i], pk, minoverlap = 1L)
    }, logical(1))
    cat(sprintf("  Reproducible-peak support: %d/%d Emerged calls confirmed; %d dropped (no strict-tier support)\n",
                sum(supported), n_emerged, sum(!supported)))
  }

  excluded_locus <- overlapsAny(gr_em, emerged_exclude_loci, minoverlap = 1L)
  if (any(excluded_locus & supported))
    cat(sprintf("  [EXCLUDE-LOCUS] %d Emerged loci overlap the construct/driver-locus exclusion region\n",
                sum(excluded_locus & supported)))
  supported <- supported & !excluded_locus

  gr_em_keep <- gr_em[supported]
  claimed_direction_keep <- claimed_direction[supported]
  n_keep <- length(gr_em_keep)

  if (n_keep > 0) {
    emerged_df <- data.frame(
      chr = as.character(seqnames(gr_em_keep)), start = start(gr_em_keep) - 1L, end = end(gr_em_keep),
      name = mcols(gr_em_keep)$name, lfc_nc14b = NA_real_, padj_nc14b = NA_real_,
      direction_nc14b = "Not differential at nc14b",
      lfc_nc14late = mcols(gr_em_keep)$log2FC, padj_nc14late = mcols(gr_em_keep)$padj,
      matched_nc14late = TRUE, fate = "Emerged",
      direction_nc14late = ifelse(claimed_direction_keep == grpA, lbl_grpA, lbl_grpB),
      stringsAsFactors = FALSE)
    fate_df <- bind_rows(fate_df, emerged_df)
  }
  cat(sprintf("  Strict Emerged calls retained: %d\n", n_keep))
}
fate_df <- as.data.frame(ungroup(fate_df))

# fate_split: same derivation as v16 (Maintained/Reversed split by which
# genotype the peak ends up open in; Converged/Emerged pass through
# unchanged -- 04_peak_fate_alluvial_plots_strict.r derives the Emerged direction split
# itself from direction_nc14late, same as it does for v16's output).
fate_df$fate_split <- dplyr::case_when(
  fate_df$fate == "Maintained" & fate_df$direction_nc14b == lbl_grpA ~ lbl_mnt_A,
  fate_df$fate == "Maintained" & fate_df$direction_nc14b == lbl_grpB ~ lbl_mnt_B,
  fate_df$fate == "Reversed"   & fate_df$direction_nc14b == lbl_grpA ~ lbl_rev_A,
  fate_df$fate == "Reversed"   & fate_df$direction_nc14b == lbl_grpB ~ lbl_rev_B,
  TRUE ~ fate_df$fate)

# resolved_breakdown: strict mode does not attempt v16's A4 Converged-
# mechanism annotation (that logic explains WHY a difference disappeared
# using the temporal strict+nearby overlap tiers, which is a separate axis
# from the significance-rescue tiers this script is designed to switch
# off). Left NA -- 04_peak_fate_alluvial_plots_strict.r detects this and skips the one
# plot that needs it, with a note, rather than failing.
fate_df$resolved_breakdown <- NA_character_

cat("\n  FINAL strict fate table:\n")
print(table(fate_df$direction_nc14b, fate_df$fate))

# Written as "peak_fate_data.csv" (not "peak_fate_data_STRICT.csv") so that
# 04_peak_fate_alluvial_plots_strict.r -- which only differs from 03_peak_fate_alluvial_plots.r
# in pointing out_contrast at this STRICT/ folder -- can read it with zero
# further changes. The STRICT/ path itself is what disambiguates this from
# the default pipeline's own peak_fate_data.csv one level up.
fn_strict_csv <- file.path(out_strict, "peak_fate_data.csv")
write.csv(fate_df, fn_strict_csv, row.names = FALSE)
cat("\nSaved:", fn_strict_csv, "\n")

strict_counts <- fate_df %>% count(fate, name = "n_strict")

################################################################################
# PART 2 -- ATTRIBUTION AGAINST v16's DEFAULT OUTPUT (if available)
################################################################################

cat("\n", strrep("=", 70), "\nPART 2: default-vs-strict comparison and rescue-tier attribution\n", strrep("=", 70), "\n\n", sep = "")

fn_default_csv <- file.path(out_contrast, "peak_fate_data.csv")

if (!file.exists(fn_default_csv)) {

  cat("  [NOTE] ", fn_default_csv, " not found -- run 01_peak_fate_classification.r first to\n",
      "  get the default-vs-strict comparison and rescue-tier attribution below.\n",
      "  Strict-only output (Part 1) has still been generated.\n", sep = "")

} else {

  default_df <- read.csv(fn_default_csv, stringsAsFactors = FALSE)
  default_counts <- default_df %>% count(fate, name = "n_default")

  comparison <- full_join(default_counts, strict_counts, by = "fate") %>%
    mutate(n_default = tidyr::replace_na(n_default, 0),
           n_strict  = tidyr::replace_na(n_strict, 0),
           n_rescued = n_default - n_strict,
           pct_rescued = ifelse(n_default > 0, round(100 * n_rescued / n_default, 1), NA)) %>%
    arrange(match(fate, c("Maintained", "Reversed", "Emerged", "Converged")))

  cat("  Default (v16, all rescue tiers) vs strict-only counts, by fate:\n\n")
  print(comparison, row.names = FALSE)

  fn_comp_csv <- file.path(out_strict, "fate_comparison_default_vs_strict.csv")
  write.csv(comparison, fn_comp_csv, row.names = FALSE)
  cat("\nSaved:", fn_comp_csv, "\n")

  # ── Attribution: which tier does each DEFAULT-run locus actually rest on? ──
  # match_confidence (Maintained/Reversed) and emerged_direction_support /
  # validation_source (Emerged) already record this per-locus in v16's
  # output -- no recomputation needed, just tabulation.
  attribution_lines <- c(
    "REPORTING NOTE FOR METHODS/SUPPLEMENT",
    strrep("=", 60), "",
    sprintf("Default pipeline (all tiers): %d classified loci (Maintained+Reversed+Converged+Emerged)",
            nrow(default_df)),
    sprintf("Strict-only (significant DAR / reproducible peak, both timepoints): %d classified loci",
            nrow(fate_df)),
    ""
  )

  if ("match_confidence" %in% colnames(default_df)) {
    mc_tab <- default_df %>%
      filter(fate %in% c("Maintained", "Reversed")) %>%
      count(match_confidence, sort = TRUE) %>%
      mutate(pct = round(100 * n / sum(n), 1))
    attribution_lines <- c(attribution_lines,
      "Maintained + Reversed calls, by evidence tier (match_confidence):", "")
    for (i in seq_len(nrow(mc_tab)))
      attribution_lines <- c(attribution_lines,
        sprintf("  %-45s %6d  (%.1f%%)", mc_tab$match_confidence[i], mc_tab$n[i], mc_tab$pct[i]))
    attribution_lines <- c(attribution_lines, "")
  }

  if ("emerged_direction_support" %in% colnames(default_df)) {
    es_tab <- default_df %>%
      filter(fate == "Emerged") %>%
      mutate(support_tier = case_when(
        grepl("^supported by reproducible", emerged_direction_support) ~ "reproducible called peak",
        grepl("^supported by called peak \\(union", emerged_direction_support) ~ "union (single-rep) rescue",
        grepl("^supported by depth", emerged_direction_support) ~ "depth+magnitude rescue",
        grepl("^UNSUPPORTED", emerged_direction_support) ~ "unsupported (kept, exclude_unsupported_emerged=FALSE)",
        grepl("^not checked", emerged_direction_support) ~ "not checked (no called-peak file)",
        TRUE ~ "other")) %>%
      count(support_tier, sort = TRUE) %>%
      mutate(pct = round(100 * n / sum(n), 1))
    attribution_lines <- c(attribution_lines,
      "Emerged calls, by support tier:", "")
    for (i in seq_len(nrow(es_tab)))
      attribution_lines <- c(attribution_lines,
        sprintf("  %-55s %6d  (%.1f%%)", es_tab$support_tier[i], es_tab$n[i], es_tab$pct[i]))
    attribution_lines <- c(attribution_lines, "")
  }

  if ("validation_source" %in% colnames(default_df)) {
    n_hcr <- sum(!is.na(default_df$validation_source) & default_df$validation_source != "")
    if (n_hcr > 0)
      attribution_lines <- c(attribution_lines,
        sprintf("HCR-validated loci added below standard ATAC threshold: %d", n_hcr), "")
  }

  fn_attr_txt <- file.path(out_strict, "fate_comparison_default_vs_strict.txt")
  writeLines(attribution_lines, fn_attr_txt)
  cat("\n", paste(attribution_lines, collapse = "\n"), "\n", sep = "")
  cat("\nSaved:", fn_attr_txt, "\n")

  # ── Comparison bar plot ─────────────────────────────────────────────────
  plot_df <- comparison %>%
    tidyr::pivot_longer(c(n_default, n_strict), names_to = "version", values_to = "n") %>%
    mutate(version = recode(version, n_default = "Default (all tiers)", n_strict = "Strict (DAR-only)"),
           fate = factor(fate, levels = c("Maintained", "Reversed", "Emerged", "Converged")))

  p <- ggplot(plot_df, aes(x = fate, y = n, fill = version)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = scales::comma(n)), position = position_dodge(width = 0.7),
              vjust = -0.4, size = 3) +
    scale_fill_manual(values = c("Default (all tiers)" = "#2c7fb8", "Strict (DAR-only)" = "#969696"),
                       name = NULL) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.12))) +
    labs(title = "Peak fate counts: default pipeline vs strict DAR-only limit",
         subtitle = "Strict = significant DAR (or reproducible called peak) required at both timepoints;\nno trend/imputation/union/depth/HCR rescue tiers",
         x = NULL, y = "Number of loci") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 8.5, colour = "grey40"))

  fn_plot <- file.path(out_strict, "fate_comparison_barplot.pdf")
  ggsave(fn_plot, p, width = 7, height = 5)
  cat("Saved:", fn_plot, "\n")
}

cat("\n", strrep("=", 70), "\nDONE. Strict-mode outputs in: ", out_strict, "\n", strrep("=", 70), "\n\n", sep = "")
