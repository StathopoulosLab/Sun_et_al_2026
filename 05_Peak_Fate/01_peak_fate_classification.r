################################################################################
# 01_peak_fate_classification.r  —  COMPREHENSIVE PEAK FATE ANALYSIS
#
# v16 CHANGE (from v15): two fixes to the per-fate gene-list CSV export.
#
# (1) Emerged loci now split into BOTv-open / BOTCv-open files, same as every
# other fate. fate_groups' geno_token assignment only ever checked
# direction_nc14b to decide BOTv vs BOTCv vs Emerged -- but Emerged loci by
# definition always have direction_nc14b=="Not differential at nc14b" (no
# nc14b DAR to have had a direction), so that check could never match for
# them and every Emerged locus landed in one undifferentiated "Emerged_
# Emerged" group regardless of which genotype it actually opened in at
# nc14late. Added an explicit fallback checking direction_nc14late (already
# set correctly per-locus in A5's emerged_df construction) for fate=="Emerged"
# rows specifically -- Converged/Maintained/Reversed token assignment is
# unaffected. Now produces genes_BOTv_Emerged.csv / genes_BOTCv_Emerged.csv.
#
# (2) save_gene_list() now (a) drops columns that are 100% NA/blank for a
# given export rather than writing them into every file regardless of
# relevance -- resolved_breakdown/mechanism_confidence are only ever
# populated for Converged loci, emerged_prior_nc14b_lfc/emerged_confidence/
# emerged_direction_support only for Emerged loci, so e.g.
# genes_BOTv_Reversed.csv previously carried all four as 100% NA. This was
# never a data bug (those columns are genuinely inapplicable outside their
# relevant fate types) but cluttered every other file. Core identifying
# columns (SYMBOL/geneId/coordinates/log2FC/nearby_genes) are exempt from
# dropping even if incidentally all-NA. And (b) nearby_genes' NA -- which
# means something different: "genuinely no gene within the window", real
# per-row information from annotate_nearby_genes() -- is now written out
# explicitly as "none within Nkb" instead of blank, so it can't be confused
# with the context-irrelevant NAs being dropped in (a).
#
# v15 CHANGE (from v14): run_chipseeker() now also calls a new
# annotate_nearby_genes() helper, adding a "nearby_genes" column alongside
# ChIPseeker's existing single-nearest-gene call. Empirical case:
# chr2R:17,760,437-17,760,673's nearest gene is ACOX1 (267bp), which is
# correct, but the locus is also ~1.8kb downstream of HLH54F -- a plain
# gene-symbol search for "HLH54F" over the exported gene-list CSVs silently
# missed this locus entirely, since ChIPseeker's SYMBOL column only ever
# holds the single nearest gene. nearby_genes lists every gene within
# gene_window_bp (default 10kb, passed through run_chipseeker(), independent
# of temporal_proximity_bp -- gene relevance and ATAC-signal proximity are
# different questions with different appropriate distance scales) as
# "SYMBOL(distance_bp)", closest first, e.g. "ACOX1(267);HLH54F(1821)".
# Purely additive: existing SYMBOL/annotation/distanceToTSS columns and every
# plot/threshold that reads them are untouched. Failure inside
# annotate_nearby_genes() is caught locally and only nulls out nearby_genes
# for that call, not the whole ChIPseeker annotation, so a problem here can't
# silently break annotation output that worked in v14.
#
# v14 CHANGE (from v13): A4's temporal_grpA/temporal_grpB flags, and A5's
# equivalent check for Emerged loci, used strict overlapsAny(minoverlap=1L)
# between a genotype-comparison DAR and the per-genotype temporal DAR set --
# a temporal DAR sitting even 1bp outside the DAR window counted as no
# evidence at all. Empirical case: BOTv_vs_BOTCv_nc14late_792
# (chr2R:17,760,437-17,760,673) sits 390bp from BOTCv_temporal_4186
# (chr2R:17,761,063-17,761,461) -- clearly the same accessible domain in the
# bigwig track, but temporal_grpB came back FALSE ("Neither changes over
# time") purely because the windows don't touch. Added a second, clearly-
# separate "nearby" tier (within temporal_proximity_bp, default 500bp, via
# findOverlaps' maxgap) alongside the existing strict-overlap tier -- nothing
# is loosened or removed, the strict flag means exactly what it meant before;
# the new tier and its recorded gap distance (temporal_grpA_gap_bp/
# temporal_grpB_gap_bp) just make previously-invisible near-adjacent evidence
# visible instead of silently discarding it. Also added an optional maxgap
# parameter (default 0, i.e. unchanged) to any_overlap_join_dir() and
# any_overlap_join_mean() for the same reason, in case you want A4b/A3c's
# trajectory-imputation rescue tier to use proximity matching too -- not
# invoked with a nonzero gap anywhere by default in this version.
#
# v13 CHANGE (from v12): read_full_results() now filters *_FULL_results.txt
# down to a single method (default "limma") and window-mode (default
# "NARROW") before returning, instead of returning every method/mode row
# undifferentiated. Those files carry multiple rows per genomic window (one
# per method x window-size combo), and any_overlap_join_mean() -- used by
# A3b's nc14late-trend rescue, A3c's "neither genotype shows a temporal shift" reclassification,
# and the Emerged sub-threshold-precursor check -- was silently averaging
# across all of them. That let noisy/less-trusted estimates (e.g. a broader
# EXTENDED window) cancel out a real, method-consistent signal from the
# same method that produced the significance-filtered DAR calls in the
# first place. See the comment on read_full_results() below for the
# empirical case (hkb locus) that surfaced this. Everything else in the
# script is unchanged from v12.
#
# SECTIONS:
#   PART A — Vent: BOTv vs BOTCv  nc14b → nc14late
#             ├─ Alluvial, stacked bar, temporal annotation plot
#             ├─ LFC scatter, magnitude shift, peak widths, emerged directions
#             ├─ ChIPseeker genomic annotation per fate group
#             ├─ Fisher's exact test: fate enrichment vs genomic background
#             ├─ LFC distribution tests (Wilcoxon) per fate
#             └─ Per-fate gene lists + GO-ready gene CSV exports
#
#   PART B — Non-vent: per-genotype temporal fate
#             ├─ BOT nc14b→nc14late, BOTC nc14b→nc14late, BOTC_oR nc14b→late
#             ├─ Each temporal DAR: Gaining / Losing + anchor overlap flag
#             ├─ ChIPseeker annotation per (geno × direction × overlap)
#             ├─ Fisher's test: Gaining vs Losing genomic feature enrichment
#             └─ Per-group gene lists + GO CSV exports
#
# OUTPUTS:
#   ./Overview_Plots/peak_fate/BOTv_vs_BOTCv/
#     peak_fate_alluvial.pdf
#     peak_fate_stackedbar.pdf
#     peak_fate_temporal_annot.pdf
#     peak_fate_lfc_scatter.pdf
#     peak_fate_magnitude_shift.pdf
#     peak_fate_widths.pdf
#     peak_fate_emerged_direction.pdf
#     peak_fate_chipseeker.pdf          ← NEW: annotation pies/bars per fate
#     peak_fate_stats_summary.txt       ← NEW: Fisher + Wilcoxon results
#     peak_fate_data.csv
#     gene_lists/                       ← NEW: per-fate gene CSVs for GO
#       genes_<fate>_<direction>.csv
#
#   ./Overview_Plots/peak_fate/NonVent/
#     nonvent_temporal_fate.pdf
#     nonvent_chipseeker.pdf            ← NEW
#     nonvent_stats_summary.txt         ← NEW
#     gene_lists/
#       genes_<geno>_<direction>.csv
#
# DEPENDENCIES:
#   Required  : GenomicRanges, ggplot2, dplyr, tidyr, patchwork, scales
#   ChIPseeker: ChIPseeker, TxDb.Dmelanogaster.UCSC.dm6.ensGene, org.Dm.eg.db
#   Optional  : ggalluvial  (install.packages("ggalluvial"))
################################################################################

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

has_alluvial   <- requireNamespace("ggalluvial",  quietly = TRUE)
has_chipseeker <- requireNamespace("ChIPseeker",  quietly = TRUE) &&
                  requireNamespace("TxDb.Dmelanogaster.UCSC.dm6.ensGene", quietly = TRUE) &&
                  requireNamespace("org.Dm.eg.db", quietly = TRUE)
has_excel      <- requireNamespace("openxlsx", quietly = TRUE)

if (!has_chipseeker)
  message("[NOTE] ChIPseeker / TxDb / org.Dm.eg.db not found — annotation plots skipped.\n",
          "  Install: BiocManager::install(c('ChIPseeker',\n",
          "    'TxDb.Dmelanogaster.UCSC.dm6.ensGene','org.Dm.eg.db'))")
if (!has_excel)
  message("[NOTE] openxlsx not found — gene-list CSVs will still be written individually, but\n",
          "  the combined multi-sheet Excel workbook will be skipped.\n",
          "  Install: install.packages('openxlsx')")

################################################################################
# CONFIG
################################################################################

dar_dir  <- "./ATAC_DARS"
out_dir  <- "./Overview_Plots"

# ── Vent comparison ────────────────────────────────────────────────────────────
stem_nc14b         <- "BOTv_vs_BOTCv"
stem_nc14late      <- "BOTv_vs_BOTCv_nc14late"
stem_temporal_grpA <- "BOTv_temporal"
stem_temporal_grpB <- "BOTCv_temporal"
grpA <- "BOTv";  col_grpA <- "#2ca02c"
grpB <- "BOTCv"; col_grpB <- "#e377c2"
lbl_grpA <- "BOT.v-biased"     # nc14b direction label (PI-requested)
lbl_grpB <- "BOTC.v-biased"    # nc14b direction label (PI-requested)

# ── Called-peak support files (OPTIONAL but strongly recommended) ─────────────
# "Emerged" direction (direction_nc14late) is currently assigned purely from
# the SIGN of DESeq2's log2FoldChange at nc14late -- see A5 below. For a
# genuine peak this is fine, but Emerged loci are, by definition, loci with
# no significant nc14b DAR to anchor them, which makes them the fate category
# most exposed to a specific failure mode: at very low read counts in one
# genotype, log2FC estimates get unstable and can land on either sign from
# noise alone, even when NEITHER genotype actually has a called peak there.
# The result is a locus labeled e.g. "Emerged (opens in BOTv)" that has no
# real BOTv peak at nc14b OR nc14late -- confirmed to happen empirically
# (see the tll/HLH54F-style diagnostics from past sessions).
#
# If you provide per-genotype, per-timepoint called-peak files below (MACS2
# narrowPeak, or any merged/consensus peak BED -- 3-column chr/start/end is
# enough), A5 will cross-check every Emerged call against them and flag any
# locus where the CLAIMED genotype has no called peak overlapping it at
# nc14late. This mirrors what check_hlh54f_emerged.r's Step 3 already did
# by hand for one locus; wiring it in here makes that check automatic and
# systematic instead of one-off.
#
# TODO: fill in the real paths/pattern for your called-peak files. Leave
# called_peak_dir as NULL (or leave files missing) to skip this check
# entirely -- everything below degrades gracefully, same as has_chipseeker/
# has_excel/has_alluvial, and Emerged calls will be produced exactly as
# before (unflagged).
called_peak_dir <- "./Merged_Peak_BEDs"
called_peak_file <- function(geno, timepoint, tier="reproducible") {
  # Two tiers, both real evidence, at different confidence levels:
  #   "reproducible" -- peak called independently in multiple reps (default,
  #                     highest confidence -- this is what guards against the
  #                     original near-zero-count sign-noise bug)
  #   "union"        -- peak called in AT LEAST ONE rep. Weaker than
  #                     reproducible, but categorically different from "no
  #                     peak ever called anywhere" -- a single-rep call is
  #                     still real MACS2 evidence, just not yet reproducibility-
  #                     confirmed. Used as a SECOND, separately-labeled rescue
  #                     tier below (see A5), never silently merged with
  #                     "reproducible" support.
  if (is.null(called_peak_dir)) return(NA_character_)
  suffix <- if (tier=="union") "_merged_union.bed" else "_merged_reproducible.bed"
  file.path(called_peak_dir, sprintf("%s_%s%s", geno, timepoint, suffix))
}

# When TRUE, Emerged loci flagged as unsupported (no called peak in the
# claimed genotype at nc14late) are dropped from fate_df entirely before the
# CSV is written -- since plot_peak_fate_v11_addon.r (Part C) and
# plot_alluvial_fixed_1_.r both read peak_fate_data.csv rather than
# recomputing fate themselves, dropping here is the single choke point and
# nothing else needs to change. When FALSE, unsupported loci are kept but
# annotated (new emerged_direction_support column) so you can inspect/filter
# them yourself -- e.g. useful the first time you run this, before you trust
# it enough to auto-drop.
exclude_unsupported_emerged <- TRUE

# ── Magnitude+depth rescue for Emerged loci without a called peak ─────────────
# The called-peak check above targets one specific failure mode: DESeq2's
# log2FC sign becomes essentially random noise at near-zero read counts, so
# "opens in BOTv" can be pure artifact when neither genotype has any real
# signal there. But requiring a MACS2-called, reproducible peak is a much
# blunter net than that -- it also excludes loci with a real, consistent,
# decent-depth effect that simply hasn't cleared MACS2's reproducibility bar
# (e.g. broad/diffuse domains). A minimum read-depth floor is the standard,
# principled guard against the sign-noise problem specifically (this is what
# DESeq2's own documentation recommends against near-zero baseMean), so it
# lets genuine low-power-but-real effects through without reopening the
# original bug. A locus is treated as "supported" if EITHER the called-peak
# check passes OR it clears both floors below.
emerged_min_depth      <- 10    # baseMean/AveExpr floor -- unchanged from the
                                 # Bourgon et al. 2010 (PNAS) independent-
                                 # filtering standard DESeq2 itself is built
                                 # on; this is already about as defensible a
                                 # floor as exists in the literature, so it's
                                 # left alone rather than loosened further
emerged_min_abs_lfc    <- 0.2   # loosened from 0.3 -- still a real, non-
                                 # trivial effect size (~15% fold-change),
                                 # a commonly used minimal bar in ATAC-seq DAR
                                 # work, not a threshold chosen to pass noise

# ── Explicit exclusion list: construct/driver-locus artifacts ────────────────
# Some loci are real, well-powered, high-magnitude signals that are STILL not
# biologically meaningful -- specifically, the promoter/regulatory region of
# the GAL4 driver transgene itself (e.g. HLH54F for BOTCv) will show large,
# genuine accessibility differences that reflect the CONSTRUCT's own
# expression, not endogenous DV patterning. This is categorically different
# from the near-zero-count sign-noise problem the depth/magnitude rescue
# guards against -- a construct artifact can pass every statistical bar
# cleanly, so no threshold change can filter it out. The correct fix is an
# explicit exclusion, checked independently of (and after) all three support
# tiers, so it can't be silently reopened by a future threshold adjustment.
#
# Add entries as you identify them; each is excluded regardless of which
# support tier (reproducible/union/depth+magnitude) it would otherwise pass.
emerged_exclude_loci <- GRanges(
  seqnames = c("chr2R"),
  ranges   = IRanges(start = c(17756000), end = c(17758500)),
  reason   = c("HLH54F driver transgene promoter (FBgn0022740) -- BOTCv driver locus, not endogenous DV signal")
)
# TODO: if you know the exact transgene insertion coordinates (rather than
# just the endogenous HLH54F gene model), widen/replace this entry to match --
# insertion-site accessibility can extend beyond the annotated gene body.

# ── Parameters ────────────────────────────────────────────────────────────────
min_recip_frac  <- 0.25

# [v14] A temporal DAR within this many bp of a genotype-comparison DAR
# counts as "nearby" evidence (A4/A5) -- separate from, and never overriding,
# the strict overlapsAny(minoverlap=1L) "overlapping" tier. Smaller than the
# 10kb window used for gene-proximity annotation, since this links two
# directly-comparable ATAC signals rather than inferring gene relevance;
# 500bp is roughly one to a few peak-widths in this dataset. Adjust to match
# your typical peak/window size -- this doesn't touch upstream DAR calling,
# so it's cheap to re-run with a different value.
temporal_proximity_bp <- 500

# RETIRED: "Deepened" and "Maintained" are now merged into a single
# "Maintained" fate category (see A3's classification), so this threshold no
# longer gates anything. Left defined (unused) only so pct_increase -- still
# computed per-locus in A3 as a continuous value, and still shown in the A10
# magnitude-shift plot -- has a documented origin. Safe to delete outright if
# you'd rather not carry dead config.
deepened_pct_thresh <- 0.05

# ── Non-vent temporal pairs ───────────────────────────────────────────────────
nv_fate_pairs <- list(
  list(stem        = "BOT_nc14d_vs_BOT_nc14b",
       geno        = "BOT",     tpA = "nc14b", tpB = "nc14late",
       col_gain    = "#145380", col_lose = "#85C1E9",
       anchor_stem = "BOT_vs_BOTR"),
  list(stem        = "BOTC_vs_BOTC_nc14late",
       geno        = "BOTC",    tpA = "nc14b", tpB = "nc14late",
       col_gain    = "#d45e00", col_lose = "#FAD7A0",
       anchor_stem = "BOT_vs_BOTC"),
  list(stem        = "BOTC_oR_late_vs_BOTC_oR_nc14b",
       geno        = "BOTC_oR", tpA = "nc14b", tpB = "nc14late",
       col_gain    = "#8c8c00", col_lose = "#E8F5E9",
       anchor_stem = "BOTC_oR_vs_BOTC")
)

# ── Palettes ──────────────────────────────────────────────────────────────────
fate_colors <- c(
  "Maintained" = "#c9579a",   # BOTCv-dominant -- using the former "Deepened" shade
                              # now that the two categories are merged
  "Reversed"   = "#8B0057",   # dark magenta (BOTv-dominant)
  "Converged"   = "#AAAAAA",
  "Emerged"    = "#F39C12"
)
fate_order <- c("Maintained","Reversed","Converged","Emerged")

dir_colors <- c(
  setNames(col_grpA, lbl_grpA),
  setNames(col_grpB, lbl_grpB),
  "Not differential at nc14b" = "#F39C12"
)

ann_feature_colors <- c(
  "Promoter"          = "#E74C3C",
  "5' UTR"            = "#F39C12",
  "3' UTR"            = "#F1C40F",
  "1st Exon"          = "#2ECC71",
  "Other Exon"        = "#27AE60",
  "1st Intron"        = "#3498DB",
  "Other Intron"      = "#2980B9",
  "Downstream"        = "#9B59B6",
  "Distal Intergenic" = "#BDC3C7"
)

################################################################################
# HELPERS
################################################################################

read_dar_full <- function(stem, search_dirs = dar_dir) {
  for (d in search_dirs) {
    bed_path <- file.path(d, paste0(stem, "_DARs.bed"))
    if (!file.exists(bed_path)) next
    bed <- tryCatch(read.table(bed_path, sep="\t", header=FALSE,
                               stringsAsFactors=FALSE), error=function(e) NULL)
    if (is.null(bed) || nrow(bed)==0) next
    n_bed <- nrow(bed)
    lfc  <- rep(NA_real_, n_bed)
    padj <- rep(NA_real_, n_bed)
    depth <- rep(NA_real_, n_bed)
    join_method <- "none"
    n_matched   <- 0L

    for (ext in c("_DARs_annotated.txt","_DARs_annotated_peaks.csv")) {
      ann_path <- file.path(d, paste0(stem, ext))
      if (!file.exists(ann_path)) next
      sep_ch <- if (grepl("\t", readLines(ann_path,1,warn=FALSE), fixed=TRUE)) "\t" else ","
      ann <- tryCatch(read.table(ann_path, header=TRUE, sep=sep_ch,
                                 stringsAsFactors=FALSE, quote="\"",
                                 fill=TRUE, comment.char="", row.names=NULL),
                      error=function(e) NULL)
      if (is.null(ann)) next
      colnames(ann)[1] <- gsub("^[^[:alnum:]]+","",colnames(ann)[1])
      lc <- intersect(c("log2FoldChange","log2FC","LFC"), colnames(ann))[1]
      pc <- intersect(c("padj","adj.P.Val","FDR","p_adj"),  colnames(ann))[1]
      # Depth column -- your annotated files carry this under a header that
      # looks like a concatenation artifact ("baseMeaAveExpr", no separating
      # tab between what should be two columns: DESeq2's baseMean and
      # limma's AveExpr). Whatever the header origin, the underlying value is
      # a real per-locus depth/expression-level signal, used below to guard
      # the Emerged magnitude+depth rescue against the exact near-zero-count
      # sign-noise failure mode the called-peak check was built to catch.
      dc <- intersect(c("baseMeaAveExpr","baseMean","AveExpr","baseMean.AveExpr"), colnames(ann))[1]
      if (is.na(lc)) next

      bed_name <- if (ncol(bed)>=4) as.character(bed[[4]]) else NULL

      # ── Tier 1: match by peak name/ID column, if one exists in both ──────
      id_col <- intersect(c("name","peak_name","peakID","peak_id","peak","ID","Row.names"),
                          colnames(ann))[1]
      if (!is.na(id_col) && !is.null(bed_name)) {
        idx <- match(bed_name, as.character(ann[[id_col]]))
        hit <- !is.na(idx)
        if (sum(hit) >= 0.5*n_bed) {   # only trust this if it actually covers most rows
          lfc[hit]  <- as.numeric(ann[[lc]])[idx[hit]]
          if (!is.na(pc)) padj[hit] <- as.numeric(ann[[pc]])[idx[hit]]
          if (!is.na(dc)) depth[hit] <- as.numeric(ann[[dc]])[idx[hit]]
          join_method <- paste0("name column '", id_col, "'")
          n_matched   <- sum(hit)
          break
        }

        # ── Tier 1b: same name column, but normalized ─────────────────────
        # Some contrasts write peak IDs as "<stem>_DAR_<n>" in the annotated
        # file while the .bed file that same contrast uses "<stem>_<n>" (no
        # "_DAR" token) -- e.g. BOTv_vs_BOTCv_DAR_1 vs BOTv_vs_BOTCv_1. This
        # is inconsistent even across contrasts (nc14late writes both sides
        # identically), so it isn't something to special-case by contrast
        # name -- just strip the literal "_DAR" token from both sides and
        # require an exact match on what's left. Still a strict string
        # match, not a positional/ordering assumption.
        norm <- function(x) gsub("_DAR(?=_|$)", "", x, perl=TRUE)
        idx_n <- match(norm(bed_name), norm(as.character(ann[[id_col]])))
        hit_n <- !is.na(idx_n)
        if (sum(hit_n) >= 0.5*n_bed) {
          lfc[hit_n]  <- as.numeric(ann[[lc]])[idx_n[hit_n]]
          if (!is.na(pc)) padj[hit_n] <- as.numeric(ann[[pc]])[idx_n[hit_n]]
          if (!is.na(dc)) depth[hit_n] <- as.numeric(ann[[dc]])[idx_n[hit_n]]
          join_method <- paste0("name column '", id_col, "' (normalized, _DAR token stripped)")
          n_matched   <- sum(hit_n)
          break
        }
      }

      # ── Tier 2: match by genomic coordinates (chr:start-end key) ─────────
      chr_col   <- intersect(c("chr","chrom","seqnames","Chr","Chrom"), colnames(ann))[1]
      start_col <- intersect(c("start","Start","chromStart"), colnames(ann))[1]
      end_col   <- intersect(c("end","End","chromEnd"), colnames(ann))[1]
      if (!is.na(chr_col) && !is.na(start_col) && !is.na(end_col)) {
        bed_key <- paste(bed[[1]], bed[[2]], bed[[3]], sep=":")
        ann_key <- paste(ann[[chr_col]], ann[[start_col]], ann[[end_col]], sep=":")
        idx <- match(bed_key, ann_key)
        hit <- !is.na(idx)
        if (sum(hit) >= 0.5*n_bed) {
          lfc[hit]  <- as.numeric(ann[[lc]])[idx[hit]]
          if (!is.na(pc)) padj[hit] <- as.numeric(ann[[pc]])[idx[hit]]
          if (!is.na(dc)) depth[hit] <- as.numeric(ann[[dc]])[idx[hit]]
          join_method <- "genomic coordinates"
          n_matched   <- sum(hit)
          break
        }
        # coordinates might be 0- vs 1-based off by one between bed and ann;
        # retry shifting the bed start by +1 before giving up on this tier
        bed_key_shift <- paste(bed[[1]], bed[[2]]+1L, bed[[3]], sep=":")
        idx2 <- match(bed_key_shift, ann_key)
        hit2 <- !is.na(idx2)
        if (sum(hit2) >= 0.5*n_bed) {
          lfc[hit2]  <- as.numeric(ann[[lc]])[idx2[hit2]]
          if (!is.na(pc)) padj[hit2] <- as.numeric(ann[[pc]])[idx2[hit2]]
          if (!is.na(dc)) depth[hit2] <- as.numeric(ann[[dc]])[idx2[hit2]]
          join_method <- "genomic coordinates (+1 start offset)"
          n_matched   <- sum(hit2)
          break
        }
      }

      # ── Tier 3: safe positional match, ONLY if row counts actually agree ──
      if (nrow(ann) == n_bed) {
        lfc  <- as.numeric(ann[[lc]])
        if (!is.na(pc)) padj <- as.numeric(ann[[pc]])
        if (!is.na(dc)) depth <- as.numeric(ann[[dc]])
        join_method <- "row position (row counts matched exactly)"
        n_matched   <- n_bed
        break
      }

      # ── Otherwise: refuse to guess. Loud warning, not silent NAs. ─────────
      warning(sprintf(
        paste0("[read_dar_full: %s] Could not safely join %s (%d rows) to %s (%d rows) -- ",
               "no usable name/ID column, no usable chr/start/end columns, and row counts ",
               "don't match (can't assume positional alignment). log2FC/padj will be NA for ",
               "ALL %d rows of this contrast. Check whether %s was filtered/deduplicated ",
               "relative to the .bed file, and add a shared name or coordinate column so a ",
               "proper join is possible."),
        stem, basename(bed_path), n_bed, basename(ann_path), nrow(ann), n_bed, basename(ann_path)))
      join_method <- "FAILED -- all NA"
      n_matched   <- 0L
    }

    cat(sprintf("  [%s] log2FC/padj joined via %s (%d/%d rows matched)\n",
               stem, join_method, n_matched, n_bed))

    if (all(is.na(lfc)) && ncol(bed)>=5) {
      raw5 <- suppressWarnings(as.numeric(bed[,5]))
      if (any(raw5<0, na.rm=TRUE)) lfc <- raw5
    }
    out <- data.frame(
      chr=as.character(bed[[1]]), start=as.integer(bed[[2]]),
      end=as.integer(bed[[3]]),
      name=if(ncol(bed)>=4) as.character(bed[[4]]) else paste0(stem,"_",seq_len(nrow(bed))),
      log2FC=as.numeric(lfc), padj=as.numeric(padj), depth=as.numeric(depth),
      stringsAsFactors=FALSE)
    stopifnot(is.data.frame(out), !inherits(out,"tbl"))
    return(out)
  }
  message("  [missing] ", stem); NULL
}

bed_to_gr <- function(df) {
  GRanges(seqnames=df$chr,
          ranges=IRanges(start=df$start+1L, end=df$end),
          name=df$name, log2FC=df$log2FC, padj=df$padj,
          depth=if ("depth" %in% names(df)) df$depth else NA_real_)
}

best_recip_join <- function(gr_query, gr_subject, min_recip=0.25) {
  hits <- findOverlaps(gr_query, gr_subject, minoverlap=1L)
  if (length(hits)==0)
    return(data.frame(i_q=integer(), i_s=integer(), recip_frac=numeric()))
  ov_w  <- width(pintersect(gr_query[queryHits(hits)], gr_subject[subjectHits(hits)]))
  recip <- pmin(ov_w/width(gr_query[queryHits(hits)]),
                ov_w/width(gr_subject[subjectHits(hits)]))
  hdf <- data.frame(i_q=queryHits(hits), i_s=subjectHits(hits), recip_frac=recip)
  hdf <- hdf[hdf$recip_frac>=min_recip,]
  if (nrow(hdf)==0) return(hdf)
  hdf %>% group_by(i_q) %>% slice_max(recip_frac,n=1,with_ties=FALSE) %>%
    ungroup() %>% as.data.frame()
}

# Any-overlap join that returns the *value* (log2FC/padj) of the best match
# in gr_subject for each locus in gr_query, rather than just T/F presence.
# Used to pull direction out of the per-genotype temporal DAR sets
# (BOTv_temporal / BOTCv_temporal) for loci whose BOTv-vs-BOTCv difference
# was called "Converged" -- presence/absence alone (the old A4 logic) can't
# tell you whether the genotype that moved was BOTv or BOTCv, or which way.
# No reciprocal-overlap requirement here (unlike best_recip_join above):
# temporal DAR widths and geno-comparison DAR widths come from different
# contrasts and aren't expected to match boundary-for-boundary, so we take
# the most significant (lowest padj) overlapping call instead.
# [v14] maxgap parameter added, default -1L (findOverlaps' own true default
# sentinel for "off" -- NOT 0L. 0L is itself a non-default maxgap value, so
# using it as this wrapper's default silently reintroduced the same maxgap+
# minoverlap conflict this comment already warns about, just one level
# removed: at the default call (no maxgap passed by the caller), findOverlaps
# was receiving an explicit maxgap=0L (non-default) together with
# minoverlap=1L (non-default) simultaneously. Using -1L as the wrapper's
# default reproduces findOverlaps' actual out-of-the-box behavior exactly,
# matching v13.) Pass maxgap=temporal_proximity_bp (or any other non-negative
# value) at a specific call site if you want THAT call to also pick up
# nearby, not just overlapping, subject ranges -- no call site does this by
# default in v15, so behavior is unchanged unless you opt in explicitly.
# minoverlap switches to 0L whenever a non-negative maxgap is requested --
# GenomicRanges' findOverlaps refuses to have BOTH maxgap and minoverlap set
# away from their own defaults simultaneously (errors with "at least one of
# 'maxgap' and 'minoverlap' must be set to its default value"); minoverlap=1L
# stays in place for the default maxgap=-1L case, matching v13's original
# strict-overlap behavior exactly.
any_overlap_join_dir <- function(gr_query, gr_subject, maxgap=-1L) {
  mo <- if (maxgap < 0L) 1L else 0L
  hits <- findOverlaps(gr_query, gr_subject, minoverlap=mo, maxgap=maxgap)
  if (length(hits)==0)
    return(data.frame(i_q=integer(), lfc=numeric(), padj=numeric()))
  hdf <- data.frame(i_q=queryHits(hits),
                    lfc=mcols(gr_subject[subjectHits(hits)])$log2FC,
                    padj=mcols(gr_subject[subjectHits(hits)])$padj)
  hdf %>% group_by(i_q) %>% slice_min(padj, n=1, with_ties=FALSE) %>%
    ungroup() %>% as.data.frame()
}

# Any-overlap join that AVERAGES log2FC across overlapping hits instead of
# picking by lowest padj. Used for the UNFILTERED results fallback below,
# where padj is frequently NA (DESeq2 sets padj=NA for independent-filtering-
# excluded rows) and isn't a meaningful tie-breaker anyway, since we're not
# filtering by significance at all here -- just recovering raw direction.
# [v14] maxgap parameter added, same rationale and same minoverlap handling
# as any_overlap_join_dir above -- default -1L, NOT 0L.
any_overlap_join_mean <- function(gr_query, gr_subject, maxgap=-1L) {
  mo <- if (maxgap < 0L) 1L else 0L
  hits <- findOverlaps(gr_query, gr_subject, minoverlap=mo, maxgap=maxgap)
  if (length(hits)==0)
    return(data.frame(i_q=integer(), lfc=numeric()))
  hdf <- data.frame(i_q=queryHits(hits),
                    lfc=mcols(gr_subject[subjectHits(hits)])$log2FC)
  hdf %>% group_by(i_q) %>% summarise(lfc=mean(lfc, na.rm=TRUE), .groups="drop") %>%
    as.data.frame()
}

# Loader for the UNFILTERED (pre-significance-filter) DESeq2/limma results
# table for a contrast -- e.g. "BOTv_temporal_FULL_results.txt" -- as opposed
# to read_dar_full() above, which reads the significance-filtered *_DARs.bed/
# *_DARs_annotated.txt pair. Generating this file requires adding an export
# step to whatever script runs the DESeq2 contrast (see write_full_results()
# in the Step1 export snippet); this loader just reads it back in if present,
# and returns NULL (with a one-line note, not an error) if it isn't -- every
# caller below is written to degrade gracefully when this file is missing,
# same pattern as has_chipseeker/has_excel/has_alluvial.
#
# Expected columns: chr, start, end, name, log2FoldChange, pvalue, padj
# (padj may be NA for many/most rows -- that's expected and fine, direction
# comes from log2FoldChange, not from padj, for this fallback).
#
# FIX (v13): the *_FULL_results.txt files carry MULTIPLE rows per genomic
# window -- one per method/mode combination (e.g. NARROW_DESeq2,
# NARROW_limma, EXTENDED200_DESeq2, EXTENDED200_limma). Previously this
# loader returned all of them undifferentiated, and any_overlap_join_mean()
# (used by A3b/A3c/emerged-precursor-check) silently AVERAGED log2FC across
# every method/window-size that happened to overlap a query locus. That let
# noisy, less-trusted estimates (e.g. a broader EXTENDED window, or the
# method NOT used for the significance-filtered DAR calls) cancel out a
# real, method-consistent signal -- confirmed empirically at the hkb locus,
# where averaging in two EXTENDED200 rows (+0.964, +0.344) against two
# NARROW rows (-1.244, -1.129) pulled a real -1.244 limma-NARROW trend down
# to a non-significant-looking -0.266 mean, and that's what kept a real
# nc14late signal from rescuing the locus out of "Converged".
#
# Fix: filter down to ONE method + ONE window-resolution before returning,
# by default the method/mode that actually produced the significance-
# filtered DAR bed files (verified against peak_fate_data.csv: its lfc
# values match NARROW_limma exactly). Every downstream caller of this
# function now compares like with like. Override filter_method/filter_mode
# per-call if a different contrast was built with a different method.
read_full_results <- function(stem, search_dirs = dar_dir,
                              filter_method = "limma", filter_mode = "NARROW") {
  for (d in search_dirs) {
    fn <- file.path(d, paste0(stem, "_FULL_results.txt"))
    if (!file.exists(fn)) next
    df <- tryCatch(read.table(fn, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                              quote="\"", fill=TRUE, comment.char=""),
                   error=function(e) NULL)
    if (is.null(df)) { message("  [", stem, "_FULL_results.txt] failed to parse -- skipping"); next }
    required_cols <- c("chr","start","end","log2FoldChange")
    if (!all(required_cols %in% colnames(df))) {
      message("  [", stem, "_FULL_results.txt] missing expected column(s) (",
              paste(setdiff(required_cols, colnames(df)), collapse=", "),
              ") -- skipping. Check write_full_results() column names match.")
      next
    }

    n_before <- nrow(df)
    if (!is.null(filter_mode) && "mode" %in% colnames(df)) {
      df <- df[df$mode == filter_mode, ]
    } else if (!is.null(filter_mode)) {
      message("  [", stem, "_FULL_results.txt] no 'mode' column found -- cannot filter to '",
              filter_mode, "'; using all rows (risk of cross-window-size averaging downstream).")
    }
    if (!is.null(filter_method) && "method" %in% colnames(df)) {
      df <- df[grepl(filter_method, df$method, ignore.case = TRUE), ]
    } else if (!is.null(filter_method)) {
      message("  [", stem, "_FULL_results.txt] no 'method' column found -- cannot filter to '",
              filter_method, "'; using all rows (risk of cross-method averaging downstream).")
    }
    cat(sprintf("  [%s_FULL_results.txt] %d/%d rows kept after filtering to method='%s', mode='%s'\n",
               stem, nrow(df), n_before,
               if (is.null(filter_method)) "any" else filter_method,
               if (is.null(filter_mode))   "any" else filter_mode))
    if (nrow(df) == 0) {
      message("  [", stem, "_FULL_results.txt] no rows left after method/mode filter -- skipping ",
              "(check filter_method/filter_mode against this file's actual 'method'/'mode' values).")
      next
    }
    return(df)
  }
  NULL
}

# Loader for an optional called-peak file (MACS2 narrowPeak or any BED with
# chr/start/end in the first three columns) for one genotype/timepoint.
# Returns a GRanges of called peaks, or NULL if called_peak_dir is NULL, the
# file doesn't exist, or it fails to parse -- every caller degrades
# gracefully, same pattern as read_full_results().
load_called_peaks <- function(geno, timepoint, tier="reproducible") {
  fp <- called_peak_file(geno, timepoint, tier=tier)
  if (is.na(fp) || !file.exists(fp)) return(NULL)
  pk <- tryCatch(read.table(fp, sep="\t", header=FALSE, stringsAsFactors=FALSE,
                            comment.char="#"),
                error=function(e) NULL)
  if (is.null(pk) || nrow(pk)==0) return(NULL)
  tryCatch(GRanges(seqnames=pk[[1]], ranges=IRanges(pk[[2]]+1L, pk[[3]])),
           error=function(e) NULL)
}

# [v15] Window-based multi-gene annotator, supplementing ChIPseeker's
# single-nearest-gene call below. ChIPseeker reports exactly one nearest
# gene per peak -- for a locus sitting between multiple genes (e.g.
# chr2R:17,760,437-17,760,673, 267bp from ACOX1's TSS but also ~1.8kb
# downstream of HLH54F), the single nearest-gene call is correct as far as
# it goes, but means a gene-symbol search for "HLH54F" silently misses a
# locus that sits squarely in HLH54F's regulatory neighborhood. This is
# additive, not a replacement: existing SYMBOL/annotation/distanceToTSS
# columns from ChIPseeker are untouched. Adds one new column, nearby_genes,
# a semicolon-separated "SYMBOL(distance_bp)" list, closest first, e.g.
# "ACOX1(267);HLH54F(1821);CG5009(4390)".
annotate_nearby_genes <- function(gr, window_bp = 10000, txdb = NULL, org_db = NULL) {
  if (is.null(txdb)) txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
  if (is.null(org_db)) org_db <- org.Dm.eg.db

  genes_gr <- GenomicFeatures::genes(txdb)  # one row per gene, gene_id = FBgn

  fbgn <- unique(genes_gr$gene_id)
  sym_map <- tryCatch(
    AnnotationDbi::select(org_db, keys = fbgn, keytype = "ENSEMBL", columns = "SYMBOL"),
    error = function(e) NULL
  )
  if (is.null(sym_map)) {
    sym_map <- tryCatch(
      AnnotationDbi::select(org_db, keys = fbgn, keytype = "FLYBASE", columns = "SYMBOL"),
      error = function(e) data.frame(ENSEMBL = fbgn, SYMBOL = fbgn)
    )
    colnames(sym_map)[colnames(sym_map) == "FLYBASE"] <- "ENSEMBL"
  }
  sym_lookup <- setNames(sym_map$SYMBOL, sym_map$ENSEMBL)

  genes_flanked <- genes_gr
  start(genes_flanked) <- pmax(1L, start(genes_gr) - window_bp)
  end(genes_flanked)   <- end(genes_gr) + window_bp

  hits <- findOverlaps(gr, genes_flanked, ignore.strand = TRUE)
  if (length(hits) == 0) return(rep(NA_character_, length(gr)))

  q_idx <- queryHits(hits); s_idx <- subjectHits(hits)
  peak_gr <- gr[q_idx]; gene_gr <- genes_gr[s_idx]
  dist_bp <- pmax(start(gene_gr) - end(peak_gr), start(peak_gr) - end(gene_gr), 0L)

  fbgn_hit   <- gene_gr$gene_id
  symbol_hit <- ifelse(fbgn_hit %in% names(sym_lookup) & !is.na(sym_lookup[fbgn_hit]),
                        sym_lookup[fbgn_hit], fbgn_hit)

  hit_df <- data.frame(q_idx = q_idx, symbol = symbol_hit, dist_bp = dist_bp,
                        stringsAsFactors = FALSE)
  nearby_str <- hit_df %>%
    dplyr::arrange(q_idx, dist_bp) %>%
    dplyr::group_by(q_idx) %>%
    dplyr::summarise(nearby_genes = paste0(symbol, "(", dist_bp, ")", collapse = ";"),
                      .groups = "drop")

  out <- rep(NA_character_, length(gr))
  out[nearby_str$q_idx] <- nearby_str$nearby_genes
  out
}

# ChIPseeker annotation helper — returns annotated GR or NULL
run_chipseeker <- function(gr, label="peaks", gene_window_bp=10000) {
  if (!has_chipseeker) return(NULL)
  suppressPackageStartupMessages({
    library(ChIPseeker)
    library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
    library(org.Dm.eg.db)
  })
  txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
  tryCatch({
    peakAnno <- annotatePeak(gr, tssRegion=c(-1000,100),
                              TxDb=txdb, annoDb="org.Dm.eg.db",
                              verbose=FALSE)
    ann_df <- as.data.frame(peakAnno)
    # [v15] additive multi-gene window column -- failure here should not
    # break the existing nearest-gene annotation, so it's wrapped separately
    ann_df$nearby_genes <- tryCatch(
      annotate_nearby_genes(gr, window_bp=gene_window_bp, txdb=txdb, org_db=org.Dm.eg.db),
      error=function(e) {
        message("  [nearby_genes warning for ", label, "]: ", e$message)
        NA_character_
      }
    )
    ann_df
  }, error=function(e) {
    message("  [ChIPseeker warning for ", label, "]: ", e$message)
    NULL
  })
}

# Collapse annotation into broad categories
simplify_annotation <- function(ann_col) {
  dplyr::case_when(
    grepl("Promoter",    ann_col, ignore.case=TRUE) ~ "Promoter",
    grepl("5' UTR",      ann_col, fixed=TRUE)       ~ "5' UTR",
    grepl("3' UTR",      ann_col, fixed=TRUE)       ~ "3' UTR",
    grepl("1st Exon",    ann_col, ignore.case=TRUE)  ~ "1st Exon",
    grepl("Exon",        ann_col, ignore.case=TRUE)  ~ "Other Exon",
    grepl("1st Intron",  ann_col, ignore.case=TRUE)  ~ "1st Intron",
    grepl("Intron",      ann_col, ignore.case=TRUE)  ~ "Other Intron",
    grepl("Downstream",  ann_col, ignore.case=TRUE)  ~ "Downstream",
    TRUE                                             ~ "Distal Intergenic"
  )
}

# Save a gene list CSV for GO analysis
# [v16] Drop columns that are 100% NA/blank for THIS export, rather than
# writing every possible column into every file regardless of relevance.
# resolved_breakdown/mechanism_confidence are only ever populated for
# Converged loci; emerged_prior_nc14b_lfc/emerged_confidence/
# emerged_direction_support only for Emerged loci -- e.g. genes_BOTv_
# Reversed.csv previously carried all four as 537/537 NA rows, which isn't a
# data bug (those columns are genuinely inapplicable to a Reversed-fate
# export) but does clutter every non-Converged/non-Emerged file with dead
# columns. Core identifying columns (SYMBOL, geneId, coordinates, log2FC)
# are exempt from dropping even if incidentally all-NA, since downstream
# consumers may reasonably expect them to always be present.
ALWAYS_KEEP_COLS <- c("SYMBOL","geneId","GENENAME","annotation","distanceToTSS",
                      "nearby_genes","seqnames","start","end","log2FC")

save_gene_list <- function(ann_df, group_label, out_gene_dir, gene_window_bp=10000) {
  if (is.null(ann_df) || nrow(ann_df)==0) return(invisible(NULL))
  cols <- intersect(c("SYMBOL","geneId","GENENAME","annotation","distanceToTSS",
                       "nearby_genes",   # [v15]
                       "seqnames","start","end","log2FC","padj",
                       "resolved_breakdown","mechanism_confidence","match_confidence",
                       "emerged_prior_nc14b_lfc","emerged_confidence",
                       "emerged_direction_support"), colnames(ann_df))
  gene_df <- ann_df[, cols, drop=FALSE]
  gene_df <- gene_df[!is.na(gene_df$SYMBOL) & gene_df$SYMBOL != "", , drop=FALSE]
  gene_df <- dplyr::distinct(gene_df, SYMBOL, .keep_all=TRUE)

  # [v16] nearby_genes NA is a DIFFERENT kind of NA from the context-specific
  # columns below -- it means "genuinely no gene within the window", real
  # per-row information, not "this column doesn't apply to this file". Made
  # explicit rather than dropped, so it isn't mistaken for missing data.
  if ("nearby_genes" %in% colnames(gene_df)) {
    gene_df$nearby_genes[is.na(gene_df$nearby_genes)] <-
      sprintf("none within %dkb", as.integer(gene_window_bp/1000))
  }

  # [v16] drop columns that are 100% NA/blank for this specific export,
  # excluding the always-kept core identifying columns above
  droppable <- setdiff(colnames(gene_df), ALWAYS_KEEP_COLS)
  all_na <- sapply(gene_df[, droppable, drop=FALSE],
                   function(x) all(is.na(x) | x==""))
  if (any(all_na)) gene_df <- gene_df[, !(colnames(gene_df) %in% names(all_na)[all_na]), drop=FALSE]

  fn <- file.path(out_gene_dir,
                  paste0("genes_", gsub("[^A-Za-z0-9_]","_",group_label), ".csv"))
  write.csv(gene_df, fn, row.names=FALSE)
  cat(sprintf("    Gene list: %s  (%d unique genes)\n", basename(fn), nrow(gene_df)))
  invisible(gene_df)
}

# ChIPseeker annotation bar plot for multiple groups
plot_annotation_bars <- function(ann_list, title="Genomic annotation by group") {
  # ann_list: named list of data.frames each with $annotation column
  rows <- lapply(names(ann_list), function(nm) {
    df <- ann_list[[nm]]
    if (is.null(df) || nrow(df)==0) return(NULL)
    df$ann_simple <- simplify_annotation(df$annotation)
    df %>% dplyr::count(ann_simple, name="n") %>%
      mutate(group=nm, pct=100*n/sum(n))
  })
  rows <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(rows) || nrow(rows)==0) return(NULL)

  feat_lvls <- c("Promoter","5' UTR","1st Exon","Other Exon",
                 "1st Intron","Other Intron","3' UTR","Downstream","Distal Intergenic")
  rows$ann_simple <- factor(rows$ann_simple, levels=rev(feat_lvls))

  ggplot(rows, aes(x=group, y=pct, fill=ann_simple)) +
    geom_bar(stat="identity", width=0.7, colour="white", linewidth=0.3) +
    geom_text(aes(label=ifelse(pct>=4, paste0(round(pct),"%" ),"")),
              position=position_stack(vjust=0.5),
              size=2.5, colour="white", fontface="bold") +
    scale_fill_manual(values=ann_feature_colors,
                      breaks=rev(feat_lvls), name="Genomic feature") +
    scale_y_continuous(labels=percent_format(scale=1), expand=c(0,0), limits=c(0,102)) +
    scale_x_discrete(labels=function(x) gsub("_","\n",x)) +
    coord_flip() +
    labs(title=title, x=NULL, y="% of peaks") +
    theme_minimal(base_size=10) +
    theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
          plot.title=element_text(size=11,face="bold"),
          legend.text=element_text(size=8), legend.title=element_text(size=8.5,face="bold"),
          legend.key.size=unit(0.75,"lines"))
}

# Fisher's exact test: feature enrichment in group vs all other peaks
fisher_feature_test <- function(ann_df, feature_pattern="Promoter",
                                background_df, group_label) {
  if (is.null(ann_df) || is.null(background_df)) return(NULL)
  in_feat  <- sum(grepl(feature_pattern, ann_df$annotation, ignore.case=TRUE))
  in_bg    <- sum(grepl(feature_pattern, background_df$annotation, ignore.case=TRUE))
  n_grp    <- nrow(ann_df)
  n_bg     <- nrow(background_df)
  mat <- matrix(c(in_feat, n_grp-in_feat,
                  in_bg-in_feat, n_bg-n_grp-(in_bg-in_feat)), 2, 2)
  ft <- tryCatch(fisher.test(mat, alternative="greater"), error=function(e) NULL)
  if (is.null(ft)) return(NULL)
  data.frame(group=group_label, feature=feature_pattern,
             n_group=n_grp, n_feature_in_group=in_feat,
             pct_feature=round(100*in_feat/n_grp,1),
             n_background=n_bg, n_feature_in_bg=in_bg,
             pct_feature_bg=round(100*in_bg/n_bg,1),
             odds_ratio=round(ft$estimate,2),
             p_value=signif(ft$p.value,3),
             stringsAsFactors=FALSE)
}

################################################################################
# ── PART A: VENT ANALYSIS ─────────────────────────────────────────────────────
################################################################################

cat("\n", strrep("=",70), "\n", sep="")
cat("PART A: VENT — BOTv vs BOTCv  nc14b → nc14late\n")
cat(strrep("=",70), "\n\n")

out_contrast     <- file.path(out_dir, "peak_fate", "BOTv_vs_BOTCv")
out_gene_dir <- file.path(out_contrast, "gene_lists")
dir.create(out_gene_dir, recursive=TRUE, showWarnings=FALSE)

# ── A1. Load ──────────────────────────────────────────────────────────────────

df_nc14b    <- read_dar_full(stem_nc14b)
df_nc14late <- read_dar_full(stem_nc14late)
if (is.null(df_nc14b) || is.null(df_nc14late))
  stop("Cannot find required BED files for vent comparison")

cat(sprintf("  nc14b genotype DARs   : %d\n", nrow(df_nc14b)))
cat(sprintf("  nc14late genotype DARs: %d\n", nrow(df_nc14late)))

df_temp_grpA <- tryCatch(read_dar_full(stem_temporal_grpA), error=function(e) NULL)
df_temp_grpB <- tryCatch(read_dar_full(stem_temporal_grpB), error=function(e) NULL)
if (!is.null(df_temp_grpA)) cat(sprintf("  %s temporal DARs: %d\n", grpA, nrow(df_temp_grpA)))
if (!is.null(df_temp_grpB)) cat(sprintf("  %s temporal DARs: %d\n", grpB, nrow(df_temp_grpB)))

# ── A2. Direction ─────────────────────────────────────────────────────────────

# BOTv_vs_BOTCv uses negated (BOTv-numerator) convention:
#   lfc > 0  →  more open in BOTv (grpA)
#   lfc < 0  →  more open in BOTCv (grpB)
df_nc14b$direction_nc14b <- ifelse(df_nc14b$log2FC>0, lbl_grpA, lbl_grpB)
df_nc14b$direction_nc14b[is.na(df_nc14b$log2FC)] <- "Unknown"

gr_nc14b    <- bed_to_gr(df_nc14b)
gr_nc14late <- bed_to_gr(df_nc14late)

# ── A3. Fate assignment ───────────────────────────────────────────────────────

join_geno <- best_recip_join(gr_nc14b, gr_nc14late, min_recip=min_recip_frac)
cat(sprintf("\n  Matched nc14b→nc14late: %d / %d (%.1f%%)\n",
            nrow(join_geno), nrow(df_nc14b),
            100*nrow(join_geno)/nrow(df_nc14b)))

fate_df <- as.data.frame(
  df_nc14b[, c("chr","start","end","name","log2FC","padj","direction_nc14b")],
  stringsAsFactors=FALSE)
colnames(fate_df)[colnames(fate_df)=="log2FC"] <- "lfc_nc14b"
colnames(fate_df)[colnames(fate_df)=="padj"]   <- "padj_nc14b"
fate_df$lfc_nc14late     <- NA_real_
fate_df$padj_nc14late    <- NA_real_
fate_df$matched_nc14late <- FALSE
fate_df$match_confidence <- NA_character_
if (nrow(join_geno)>0) {
  fate_df$lfc_nc14late[join_geno$i_q]     <- mcols(gr_nc14late[join_geno$i_s])$log2FC
  fate_df$padj_nc14late[join_geno$i_q]    <- mcols(gr_nc14late[join_geno$i_s])$padj
  fate_df$matched_nc14late[join_geno$i_q] <- TRUE
  fate_df$match_confidence[join_geno$i_q] <- "significant DAR"
}

# ── A3b. Unfiltered-trend fallback for unmatched nc14b DARs ───────────────────
# matched_nc14late above only counts a match if the locus is ALSO a
# significant BOTv-vs-BOTCv DAR at nc14late -- but gr_nc14late (the pool it's
# matched against) has only 1,371 peaks vs 6,323 at nc14b. Even a perfect
# match rate could never exceed 1,371/6,323 (21.7%) of nc14b peaks, and the
# actual match rate is far below even that ceiling -- a nearly 5-fold drop in
# significant DAR calls between two timepoints ~30 minutes apart, which looks
# far more like a power/threshold difference between the two DESeq2 tests
# than 90%+ of real genotype differences genuinely resolving that fast. Any
# locus that's really still different at nc14late but fell just short of
# nc14late's (evidently stricter) significance bar would incorrectly default
# to "Converged" here rather than Maintained/Deepened/Reversed.
#
# If BOTv_vs_BOTCv_nc14late_FULL_results.txt exists (unfiltered results, every
# tested peak regardless of significance -- see write_full_results() in the
# Step1 export snippet), use it to give currently-unmatched nc14b DARs a
# second chance: same overlap + sign-consistency logic as the temporal
# fallback above, using a genotype-difference-appropriate threshold rather
# than the temporal one. This can only ever ADD matches, never remove one
# that already passed the strict significance-based join.
nc14late_trend_thresh <- 0.5   # a log2FC-scale threshold (not percent-based like deepened_pct_thresh) -- distinguishes
                               # "still meaningfully different" from noise-level LFC

# Consistency band for the "still basically the same" rescue path below --
# same reasoning as maintained_consistency_band in tier 3 (A3c). See that
# comment block for the full explanation of why a fixed magnitude floor
# alone structurally can't ever produce a Maintained call for small-to-
# moderate starting effect sizes.
maintained_consistency_band <- c(0.7, 1.3)   # candidate within 70-130% of |lfc_nc14b|

full_nc14late <- read_full_results("BOTv_vs_BOTCv_nc14late")
if (is.null(full_nc14late)) {
  cat("  [NOTE] BOTv_vs_BOTCv_nc14late_FULL_results.txt not found -- unmatched nc14b\n",
      "  DARs default straight to 'Converged'. See write_full_results() in the Step1\n",
      "  export snippet to let trending-but-not-significant matches through instead.\n")
} else {
  unmatched_idx <- which(!fate_df$matched_nc14late)
  if (length(unmatched_idx) > 0) {
    gr_unmatched <- GRanges(seqnames=fate_df$chr[unmatched_idx],
                            ranges=IRanges(fate_df$start[unmatched_idx]+1L,
                                          fate_df$end[unmatched_idx]))
    gr_full_nc14late <- GRanges(seqnames=full_nc14late$chr,
                                ranges=IRanges(full_nc14late$start+1L, full_nc14late$end),
                                log2FC=full_nc14late$log2FoldChange)
    j_trend <- any_overlap_join_mean(gr_unmatched, gr_full_nc14late)
    if (nrow(j_trend) > 0) {
      trend_lfc <- j_trend$lfc
      lfc_nc14b_unm <- fate_df$lfc_nc14b[unmatched_idx[j_trend$i_q]]
      # Path 1 (unchanged): candidate clears an absolute magnitude floor,
      # regardless of sign relative to lfc_nc14b -- catches real Deepened/
      # Reversed trends. NOT restricted to same sign; the downstream case_when
      # already sorts same-sign vs opposite-sign correctly on its own.
      trend_ok_magnitude <- !is.na(trend_lfc) & abs(trend_lfc) > nc14late_trend_thresh
      # Path 2 (new): candidate is close to lfc_nc14b in relative magnitude
      # AND same sign -- direct evidence of persistence/stability, valid even
      # when both values are individually below the absolute floor above.
      # Without this, tiers 2/3 can only ever rescue Deepened/Reversed calls:
      # clearing a fixed floor from a small starting lfc_nc14b is ALWAYS a
      # >5% relative jump (e.g. 0.4->0.5 is +7%), so a small/moderate-effect
      # peak that's genuinely just holding steady could never pass the
      # magnitude-only test without ALSO tripping the Deepened threshold.
      trend_ok_maintained <- !is.na(trend_lfc) &
        sign(trend_lfc)==sign(lfc_nc14b_unm) &
        abs(trend_lfc) >= maintained_consistency_band[1]*abs(lfc_nc14b_unm) &
        abs(trend_lfc) <= maintained_consistency_band[2]*abs(lfc_nc14b_unm)
      trend_ok <- trend_ok_magnitude | trend_ok_maintained
      rescued <- unmatched_idx[j_trend$i_q[trend_ok]]
      if (length(rescued) > 0) {
        fate_df$lfc_nc14late[rescued]     <- trend_lfc[trend_ok]
        fate_df$matched_nc14late[rescued] <- TRUE
        fate_df$match_confidence[rescued] <- ifelse(
          trend_ok_maintained[trend_ok] & !trend_ok_magnitude[trend_ok],
          "trending (stability-consistency rescue)",
          "trending, not FDR-significant")
      }
      cat(sprintf("  Unfiltered-trend fallback: rescued %d/%d unmatched nc14b DARs (%d via magnitude >%.2f, %d via stability-consistency [%.0f-%.0f%% of nc14b])\n",
                 length(rescued), length(unmatched_idx),
                 sum(trend_ok_magnitude[trend_ok]),
                 nc14late_trend_thresh,
                 sum(trend_ok_maintained[trend_ok] & !trend_ok_magnitude[trend_ok]),
                 100*maintained_consistency_band[1], 100*maintained_consistency_band[2]))
    }
  }
}

# ── A3c. Trajectory imputation using BOTv_temporal / BOTCv_temporal ───────────
# For loci STILL unmatched after tiers 1 (significant nc14late DAR) and 2
# (unfiltered nc14late trend), we have zero direct nc14late genotype-
# comparison evidence at all -- but we may still be able to say something,
# using data that's already loaded regardless of whether the Step1 export
# snippet has been run: the SIGNIFICANT BOTv_temporal and BOTCv_temporal DAR
# calls (1,283 and 5,909 peaks respectively -- much larger, better-powered
# pools than BOTv_vs_BOTCv_nc14late's 1,371).
#
# The algebra: if lfc_nc14b approximates log2(BOTv/BOTCv) at nc14b, and each
# genotype's own accessibility shifts by some amount between nc14b and
# nc14late (that's exactly what BOTv_temporal/BOTCv_temporal measure), then
#   predicted_lfc_nc14late = lfc_nc14b + (BOTv's own shift) - (BOTCv's own shift)
# This lets a locus be reclassified out of "Converged" even with NO nc14late
# genotype-comparison signal, purely from each genotype's independently-
# measured temporal trajectory -- which is exactly the "early vs late" data
# already sitting in BOTv_temporal_DARs.bed / BOTCv_temporal_DARs.bed.
#
# Requirements, deliberately conservative:
#   - At least ONE of the two genotypes must have a SIGNIFICANT temporal DAR
#     overlapping the locus (if neither does, there's no real evidence to
#     impute from, and it correctly stays "Converged: neither genotype shows a temporal shift").
#     A missing side is treated as 0 shift, not as evidence of stability --
#     that's an assumption, not a measurement, and is why this tier's
#     match_confidence label says "imputed", not "trending" or "significant".
#   - This assumes the two contrasts add linearly (no genotype x time
#     interaction beyond what's captured by each contrast alone). That's a
#     simplification -- worth spot-checking against a locus you've already
#     inspected visually (e.g. Piezo/aop/net) before trusting it broadly.
#   - Slightly stricter magnitude threshold than tier 2, since this estimate
#     compounds uncertainty from two separate contrasts rather than measuring
#     the quantity of interest directly.
#
# Same magnitude-floor bias as tier 2 applies here too, and for the same
# reason -- see maintained_consistency_band above -- so this tier gets the
# identical second rescue path: a predicted value close to lfc_nc14b in
# relative terms counts as evidence of persistence even when it doesn't
# individually clear impute_thresh.
impute_thresh <- 0.6

still_unmatched <- which(!fate_df$matched_nc14late)
if (length(still_unmatched) > 0 && (!is.null(df_temp_grpA) || !is.null(df_temp_grpB))) {
  gr_still <- GRanges(seqnames=fate_df$chr[still_unmatched],
                      ranges=IRanges(fate_df$start[still_unmatched]+1L,
                                    fate_df$end[still_unmatched]))
  trend_a3 <- rep(NA_real_, length(still_unmatched))
  trend_b3 <- rep(NA_real_, length(still_unmatched))
  if (!is.null(df_temp_grpA) && nrow(df_temp_grpA)>0) {
    gr_temp_a3 <- bed_to_gr(df_temp_grpA)
    j_a3 <- any_overlap_join_dir(gr_still, gr_temp_a3)
    if (nrow(j_a3)>0) trend_a3[j_a3$i_q] <- j_a3$lfc
  }
  if (!is.null(df_temp_grpB) && nrow(df_temp_grpB)>0) {
    gr_temp_b3 <- bed_to_gr(df_temp_grpB)
    j_b3 <- any_overlap_join_dir(gr_still, gr_temp_b3)
    if (nrow(j_b3)>0) trend_b3[j_b3$i_q] <- j_b3$lfc
  }
  has_evidence <- !is.na(trend_a3) | !is.na(trend_b3)
  predicted_lfc <- fate_df$lfc_nc14b[still_unmatched] +
                   dplyr::coalesce(trend_a3, 0) - dplyr::coalesce(trend_b3, 0)
  lfc_nc14b_stu <- fate_df$lfc_nc14b[still_unmatched]
  impute_ok_magnitude   <- has_evidence & abs(predicted_lfc) > impute_thresh
  impute_ok_maintained  <- has_evidence &
    sign(predicted_lfc)==sign(lfc_nc14b_stu) &
    abs(predicted_lfc) >= maintained_consistency_band[1]*abs(lfc_nc14b_stu) &
    abs(predicted_lfc) <= maintained_consistency_band[2]*abs(lfc_nc14b_stu)
  impute_ok <- impute_ok_magnitude | impute_ok_maintained
  imputed <- still_unmatched[impute_ok]
  if (length(imputed) > 0) {
    fate_df$lfc_nc14late[imputed]     <- predicted_lfc[impute_ok]
    fate_df$matched_nc14late[imputed] <- TRUE
    fate_df$match_confidence[imputed] <- ifelse(
      impute_ok_maintained[impute_ok] & !impute_ok_magnitude[impute_ok],
      "imputed (stability-consistency rescue)",
      "imputed from BOTv/BOTCv temporal trajectories")
  }
  cat(sprintf("  Trajectory-imputation fallback: rescued %d/%d still-unmatched nc14b DARs (%d via magnitude >%.2f, %d via stability-consistency [%.0f-%.0f%% of nc14b])\n",
             length(imputed), length(still_unmatched),
             sum(impute_ok_magnitude[impute_ok]), impute_thresh,
             sum(impute_ok_maintained[impute_ok] & !impute_ok_magnitude[impute_ok]),
             100*maintained_consistency_band[1], 100*maintained_consistency_band[2]))
}

fate_df$match_confidence[is.na(fate_df$match_confidence)] <- "no nc14late match"

fate_df <- fate_df %>% mutate(
  # Fold-change magnitude of the genotype-difference at each timepoint (always
  # >=1 since it's on |log2FC|), and the relative increase between them --
  # this is what "Deepened" actually tests now, not a fixed log2FC bump.
  fc_nc14b_mag     = 2^abs(lfc_nc14b),
  fc_nc14late_mag  = 2^abs(lfc_nc14late),
  pct_increase     = (fc_nc14late_mag - fc_nc14b_mag) / fc_nc14b_mag,
  fate = case_when(
    !matched_nc14late                                          ~ "Converged",
    matched_nc14late & is.na(lfc_nc14late)                   ~ "Converged",
    matched_nc14late & !is.na(lfc_nc14late) &
      sign(lfc_nc14late)==sign(lfc_nc14b)                    ~ "Maintained",
    matched_nc14late & !is.na(lfc_nc14late) &
      sign(lfc_nc14late)!=sign(lfc_nc14b)                    ~ "Reversed",
    TRUE                                                       ~ "Converged"),
  direction_nc14late = case_when(
    fate=="Converged"                           ~ "Converged",
    !is.na(lfc_nc14late) & lfc_nc14late>0    ~ lbl_grpA,
    !is.na(lfc_nc14late) & lfc_nc14late<0    ~ lbl_grpB,
    TRUE                                       ~ "Converged"))

cat("\n  Fate breakdown:\n")
print(table(fate_df$direction_nc14b, fate_df$fate))

# ── A4. Temporal annotation ───────────────────────────────────────────────────

fate_df$temporal_grpA        <- FALSE   # strict overlap -- UNCHANGED meaning from v13
fate_df$temporal_grpB        <- FALSE
fate_df$temporal_grpA_nearby <- FALSE   # [v14] within temporal_proximity_bp, not overlapping
fate_df$temporal_grpB_nearby <- FALSE
fate_df$temporal_grpA_gap_bp <- NA_integer_  # [v14] actual gap distance, for inspection
fate_df$temporal_grpB_gap_bp <- NA_integer_

if (!is.null(df_temp_grpA) && nrow(df_temp_grpA)>0) {
  gr_ta <- bed_to_gr(df_temp_grpA)
  fate_df$temporal_grpA <- overlapsAny(gr_nc14b, gr_ta, minoverlap=1L)
  # [v14] "nearby but not overlapping": within the proximity window per
  # maxgap, but the strict check above already came back FALSE
  within_window <- overlapsAny(gr_nc14b, gr_ta, maxgap=temporal_proximity_bp, minoverlap=0L)
  fate_df$temporal_grpA_nearby <- within_window & !fate_df$temporal_grpA
  if (any(fate_df$temporal_grpA_nearby)) {
    idx <- which(fate_df$temporal_grpA_nearby)
    d   <- distanceToNearest(gr_nc14b[idx], gr_ta, ignore.strand=TRUE)
    fate_df$temporal_grpA_gap_bp[idx[queryHits(d)]] <- mcols(d)$distance
  }
}
if (!is.null(df_temp_grpB) && nrow(df_temp_grpB)>0) {
  gr_tb <- bed_to_gr(df_temp_grpB)
  fate_df$temporal_grpB <- overlapsAny(gr_nc14b, gr_tb, minoverlap=1L)
  within_window <- overlapsAny(gr_nc14b, gr_tb, maxgap=temporal_proximity_bp, minoverlap=0L)
  fate_df$temporal_grpB_nearby <- within_window & !fate_df$temporal_grpB
  if (any(fate_df$temporal_grpB_nearby)) {
    idx <- which(fate_df$temporal_grpB_nearby)
    d   <- distanceToNearest(gr_nc14b[idx], gr_tb, ignore.strand=TRUE)
    fate_df$temporal_grpB_gap_bp[idx[queryHits(d)]] <- mcols(d)$distance
  }
}

# [v14] temporal_label now distinguishes overlapping vs nearby vs none, per
# genotype, instead of collapsing straight to a boolean -- collapsing that
# distinction in v13 was exactly what discarded the BOTCv_temporal_4186
# evidence for the chr2R:17,760,437 locus (390bp away, real signal, silently
# treated as "Neither changes over time").
label_tier <- function(overlap, nearby, gap, label) {
  case_when(
    overlap ~ paste0(label, " changes over time (overlapping)"),
    nearby  ~ paste0(label, " changes over time (nearby, ", gap, "bp)"),
    TRUE    ~ NA_character_
  )
}
.tier_A <- label_tier(fate_df$temporal_grpA, fate_df$temporal_grpA_nearby,
                       fate_df$temporal_grpA_gap_bp, grpA)
.tier_B <- label_tier(fate_df$temporal_grpB, fate_df$temporal_grpB_nearby,
                       fate_df$temporal_grpB_gap_bp, grpB)
fate_df$temporal_label <- case_when(
  !is.na(.tier_A) & !is.na(.tier_B) ~ paste0("Both change -- ", .tier_A, "; ", .tier_B),
  !is.na(.tier_A)                    ~ .tier_A,
  !is.na(.tier_B)                    ~ .tier_B,
  TRUE                                ~ "Neither changes over time")
rm(.tier_A, .tier_B)

# ── A4b. Resolution mechanism via per-genotype temporal DARs ──────────────────
# "Converged" (A3 above) only tells us the BOTv-vs-BOTCv difference is no
# longer a significant DAR at nc14late -- it says nothing about *why*. The
# temporal_grpA/grpB flags just above only capture presence/absence, which
# throws away exactly the information needed to explain it. Here we pull the
# actual log2FC of the best-overlapping BOTv_temporal / BOTCv_temporal DAR
# (nc14b -> nc14late within each genotype) so we can tell whether it was
# BOTv closing, BOTCv opening, both, or neither, that drove the convergence.
#
# SIGN CONVENTION -- PLEASE VERIFY: this assumes BOTv_temporal / BOTCv_temporal
# log2FC is signed nc14late-vs-nc14b (positive = more open at nc14late, i.e.
# "opens"). If your temporal contrasts were built with nc14b as the numerator
# instead, flip the ifelse() below (opens <-> closes). Sanity-check against a
# locus you've already inspected visually -- e.g. Piezo/aop/net, where BOTCv
# visibly gains signal from nc14b to nc14late -- and confirm it comes back as
# "opens" for temporal_dir_grpB before trusting the rest.

fate_df$temporal_lfc_grpA <- NA_real_
fate_df$temporal_lfc_grpB <- NA_real_
fate_df$temporal_dir_grpA <- NA_character_
fate_df$temporal_dir_grpB <- NA_character_

if (!is.null(df_temp_grpA) && nrow(df_temp_grpA)>0) {
  j_ta <- any_overlap_join_dir(gr_nc14b, gr_ta)
  if (nrow(j_ta)>0) {
    fate_df$temporal_lfc_grpA[j_ta$i_q] <- j_ta$lfc
    fate_df$temporal_dir_grpA[j_ta$i_q] <- ifelse(j_ta$lfc>0, "opens", "closes")
  }
}
if (!is.null(df_temp_grpB) && nrow(df_temp_grpB)>0) {
  j_tb <- any_overlap_join_dir(gr_nc14b, gr_tb)
  if (nrow(j_tb)>0) {
    fate_df$temporal_lfc_grpB[j_tb$i_q] <- j_tb$lfc
    fate_df$temporal_dir_grpB[j_tb$i_q] <- ifelse(j_tb$lfc>0, "opens", "closes")
  }
}

fate_df$resolution_mechanism <- NA_character_
resolved_idx <- which(fate_df$fate=="Converged")
if (length(resolved_idx)>0) {
  fate_df$resolution_mechanism[resolved_idx] <- with(fate_df[resolved_idx,], case_when(
    temporal_dir_grpA=="closes" & temporal_dir_grpB=="opens"          ~ paste0("Converged: ",grpA," closes + ",grpB," opens"),
    temporal_dir_grpA=="opens"  & temporal_dir_grpB=="closes"         ~ paste0("Converged: ",grpB," closes + ",grpA," opens"),
    is.na(temporal_dir_grpA)    & temporal_dir_grpB=="opens"          ~ paste0("Converged: ",grpB," opens (",grpA," unchanged)"),
    is.na(temporal_dir_grpA)    & temporal_dir_grpB=="closes"         ~ paste0("Converged: ",grpB," closes (",grpA," unchanged)"),
    temporal_dir_grpA=="closes" & is.na(temporal_dir_grpB)            ~ paste0("Converged: ",grpA," closes (",grpB," unchanged)"),
    temporal_dir_grpA=="opens"  & is.na(temporal_dir_grpB)            ~ paste0("Converged: ",grpA," opens (",grpB," unchanged)"),
    !is.na(temporal_dir_grpA) & temporal_dir_grpA==temporal_dir_grpB  ~ paste0("Converged: both shift same direction (",temporal_dir_grpA,")"),
    TRUE                                                               ~ "No temporal DAR detected in either genotype"))
}

cat("\n  Converged-locus mechanism breakdown (why the BOTv-vs-BOTCv difference disappeared):\n")
print(table(fate_df$resolution_mechanism, useNA="ifany"))

# ── A5. Emerged peaks ─────────────────────────────────────────────────────────

emerged_mask <- !overlapsAny(gr_nc14late, gr_nc14b, minoverlap=1L)
n_emerged    <- sum(emerged_mask)
cat(sprintf("\n  De novo peaks at nc14late: %d\n", n_emerged))

if (n_emerged>0) {
  gr_em <- gr_nc14late[emerged_mask]
  emrg_ta <- rep(FALSE, n_emerged)   # [v14] strict overlap -- UNCHANGED meaning from v13
  emrg_tb <- rep(FALSE, n_emerged)
  emrg_ta_nearby <- rep(FALSE, n_emerged)   # [v14] within temporal_proximity_bp, not overlapping
  emrg_tb_nearby <- rep(FALSE, n_emerged)
  if (!is.null(df_temp_grpA) && nrow(df_temp_grpA)>0) {
    emrg_ta <- overlapsAny(gr_em, gr_ta, minoverlap=1L)
    emrg_ta_nearby <- overlapsAny(gr_em, gr_ta, maxgap=temporal_proximity_bp, minoverlap=0L) & !emrg_ta
  }
  if (!is.null(df_temp_grpB) && nrow(df_temp_grpB)>0) {
    emrg_tb <- overlapsAny(gr_em, gr_tb, minoverlap=1L)
    emrg_tb_nearby <- overlapsAny(gr_em, gr_tb, maxgap=temporal_proximity_bp, minoverlap=0L) & !emrg_tb
  }

  # "Emerged" means no SIGNIFICANT nc14b DAR overlaps this nc14late peak --
  # but a real, just-below-threshold difference at nc14b would look
  # identical to true de novo appearance under that test alone, which
  # overstates how many peaks are genuinely appearing from nothing rather
  # than growing out of an existing weak signal. If BOTv_vs_BOTCv_FULL_
  # results.txt exists (unfiltered, every tested peak -- see
  # write_full_results() in the Step1 export snippet), check each Emerged
  # peak against the RAW nc14b value regardless of significance. This only
  # ever ANNOTATES -- it never removes a peak from "Emerged" or changes any
  # existing plot/color logic, since so much already keys off that label;
  # it just exposes the distinction for filtering/interpretation.
  emerged_prior_thresh <- 0.25   # deliberately low bar -- "any hint of prior
                                 # signal", not "was it a real difference"
  full_nc14b <- read_full_results("BOTv_vs_BOTCv")
  emerged_prior_lfc <- rep(NA_real_, n_emerged)
  if (is.null(full_nc14b)) {
    cat("  [NOTE] BOTv_vs_BOTCv_FULL_results.txt not found -- 'Emerged' peaks\n",
        "  can't be checked for a sub-threshold nc14b precursor. See\n",
        "  write_full_results() in the Step1 export snippet to enable this check.\n")
  } else {
    gr_full_nc14b <- GRanges(seqnames=full_nc14b$chr,
                             ranges=IRanges(full_nc14b$start+1L, full_nc14b$end),
                             log2FC=full_nc14b$log2FoldChange)
    j_em <- any_overlap_join_mean(gr_em, gr_full_nc14b)
    if (nrow(j_em)>0) emerged_prior_lfc[j_em$i_q] <- j_em$lfc
  }
  emerged_had_prior <- !is.na(emerged_prior_lfc) & abs(emerged_prior_lfc) > emerged_prior_thresh
  emerged_confidence <- if (is.null(full_nc14b)) {
    rep("unfiltered nc14b data unavailable", n_emerged)
  } else {
    ifelse(emerged_had_prior, "grew from sub-threshold nc14b signal", "true de novo (no nc14b evidence)")
  }
  if (!is.null(full_nc14b)) {
    cat(sprintf("  Emerged-peak check: %d/%d actually had a sub-threshold nc14b precursor (|log2FC|>%.2f) rather than being truly de novo\n",
               sum(emerged_had_prior), n_emerged, emerged_prior_thresh))
  }

  claimed_direction_nc14late <- ifelse(!is.na(mcols(gr_em)$log2FC) & mcols(gr_em)$log2FC>0,
                                       grpA, grpB)

  # ── Called-peak support check ──────────────────────────────────────────────
  # A locus labeled "Emerged (opens in BOTv)" should have a real, called BOTv
  # peak at nc14late overlapping it. If it doesn't, the positive log2FC
  # driving that label is more likely DESeq2 sign noise on a near-zero-count
  # locus than a genuine BOTv-specific signal -- variance blows up as counts
  # approach zero, and the point estimate can land on either sign essentially
  # at random. This is exactly the failure mode confirmed by hand for
  # individual loci in past sessions (check_hlh54f_emerged.r Step 3): zero
  # called peaks in the "winning" genotype at either timepoint, but a
  # positive log2FC anyway.
  #
  # Like emerged_prior_lfc/emerged_confidence above, this only ANNOTATES by
  # default (exclude_unsupported_emerged=FALSE) -- set that flag in CONFIG to
  # TRUE once you've spot-checked a few flagged loci and trust the filter.
  peaks_grpA_nc14late <- load_called_peaks(grpA, "nc14late", tier="reproducible")
  peaks_grpB_nc14late <- load_called_peaks(grpB, "nc14late", tier="reproducible")
  peaks_grpA_nc14late_union <- load_called_peaks(grpA, "nc14late", tier="union")
  peaks_grpB_nc14late_union <- load_called_peaks(grpB, "nc14late", tier="union")
  has_called_peaks <- !is.null(peaks_grpA_nc14late) || !is.null(peaks_grpB_nc14late) ||
                       !is.null(peaks_grpA_nc14late_union) || !is.null(peaks_grpB_nc14late_union)
  cat(sprintf("  Called-peak file status: %s reproducible=%s | %s reproducible=%s | %s union=%s | %s union=%s\n",
             grpA, !is.null(peaks_grpA_nc14late),
             grpB, !is.null(peaks_grpB_nc14late),
             grpA, !is.null(peaks_grpA_nc14late_union),
             grpB, !is.null(peaks_grpB_nc14late_union)))
  if (!is.null(peaks_grpA_nc14late))       cat("   ", grpA, "reproducible file:", called_peak_file(grpA,"nc14late","reproducible"), "\n")
  if (!is.null(peaks_grpA_nc14late_union)) cat("   ", grpA, "union file:       ", called_peak_file(grpA,"nc14late","union"), "\n")
  if (!is.null(peaks_grpB_nc14late))       cat("   ", grpB, "reproducible file:", called_peak_file(grpB,"nc14late","reproducible"), "\n")
  if (!is.null(peaks_grpB_nc14late_union)) cat("   ", grpB, "union file:       ", called_peak_file(grpB,"nc14late","union"), "\n")

  if (!has_called_peaks) {
    emerged_direction_support <- rep("not checked (no called-peak file configured -- see CONFIG)", n_emerged)
    cat("  [NOTE] No called-peak files found for", grpA, "/", grpB, "at nc14late --\n",
        "  Emerged direction calls are UNVERIFIED against actual called peaks.\n",
        "  See the 'Called-peak support files' CONFIG block near the top of this\n",
        "  script to enable this check.\n")
  } else {
    supported_by_peak <- logical(n_emerged)
    supported_by_union <- logical(n_emerged)
    for (i in seq_len(n_emerged)) {
      pk       <- if (claimed_direction_nc14late[i]==grpA) peaks_grpA_nc14late       else peaks_grpB_nc14late
      pk_union <- if (claimed_direction_nc14late[i]==grpA) peaks_grpA_nc14late_union else peaks_grpB_nc14late_union
      supported_by_peak[i]  <- !is.null(pk)       && overlapsAny(gr_em[i], pk,       minoverlap=1L)
      supported_by_union[i] <- !is.null(pk_union) && overlapsAny(gr_em[i], pk_union, minoverlap=1L)
    }
    # union support only counts for loci the reproducible check didn't
    # already cover -- keeps the tiers cleanly separated rather than
    # double-counting the same locus under two labels
    supported_by_union <- supported_by_union & !supported_by_peak

    # Magnitude+depth rescue -- see CONFIG block for rationale. Applied only
    # to loci NEITHER called-peak check already supported, so it can only
    # ADD coverage, never override a called-peak-based exclusion.
    depth_em <- mcols(gr_em)$depth
    lfc_em   <- mcols(gr_em)$log2FC
    supported_by_depth <- !is.na(depth_em) & depth_em >= emerged_min_depth &
                          !is.na(lfc_em)   & abs(lfc_em) >= emerged_min_abs_lfc &
                          !supported_by_peak & !supported_by_union
    if (all(is.na(depth_em))) {
      cat("  [NOTE] Depth column (baseMean/AveExpr) not available for this contrast's\n",
          "  annotated file -- magnitude+depth rescue inactive; falling back to\n",
          "  called-peak support only. See read_dar_full()'s depth-column detection\n",
          "  if you expect this to be populated.\n")
      supported_by_depth <- rep(FALSE, n_emerged)
    }

    supported <- supported_by_peak | supported_by_union | supported_by_depth

    # Explicit exclusion list overrides ALL support tiers -- a construct/
    # driver-locus artifact can be well-powered and high-magnitude, so no
    # statistical tier should be able to rescue it. See CONFIG for rationale.
    excluded_locus <- overlapsAny(gr_em, emerged_exclude_loci, minoverlap=1L)
    n_excluded_locus <- sum(excluded_locus & supported)
    if (n_excluded_locus > 0) {
      cat(sprintf("  [EXCLUDE-LOCUS] %d Emerged loci overlap an explicit exclusion region (construct/driver-locus artifact) -- excluded regardless of support tier.\n",
                 n_excluded_locus))
    }
    supported <- supported & !excluded_locus

    n_rescued_union <- sum(supported_by_union & supported)
    n_rescued_depth <- sum(supported_by_depth & supported)
    if (n_rescued_union > 0) {
      cat(sprintf("  [RESCUE-UNION] %d/%d Emerged loci lack a reproducible peak but ARE called in at least one replicate (union bed) -- kept as supported.\n",
                 n_rescued_union, n_emerged))
    }
    if (n_rescued_depth > 0) {
      cat(sprintf("  [RESCUE-DEPTH] %d/%d Emerged loci lack any called peak (reproducible or union) but clear the depth+magnitude floor (baseMean>=%.0f, |log2FC|>=%.2f) -- kept as supported.\n",
                 n_rescued_depth, n_emerged, emerged_min_depth, emerged_min_abs_lfc))
    }

    emerged_direction_support <- ifelse(
      excluded_locus, "EXCLUDED -- overlaps known construct/driver-locus artifact region (see emerged_exclude_loci in CONFIG)",
      ifelse(supported_by_peak, "supported by reproducible called peak",
      ifelse(supported_by_union, "supported by called peak (union, single-replicate)",
      ifelse(supported_by_depth,
             sprintf("supported by depth+magnitude (baseMean=%.0f, log2FC=%.2f) -- no called peak",
                     depth_em, lfc_em),
             sprintf("UNSUPPORTED -- no called %s peak (reproducible or union) at nc14late AND fails depth/magnitude floor; likely low-count instability, not a real %s-specific signal",
                     claimed_direction_nc14late, claimed_direction_nc14late)))))
    n_unsupported <- sum(!supported)
    if (n_unsupported>0) {
      cat(sprintf("  [WARNING] %d/%d Emerged calls have NO called peak (either tier) in the claimed genotype at nc14late AND fail the depth/magnitude rescue\n",
                 n_unsupported, n_emerged))
      cat("   ", if (exclude_unsupported_emerged)
            "exclude_unsupported_emerged=TRUE -- these will be dropped from fate_df below.\n"
          else
            "exclude_unsupported_emerged=FALSE -- kept, flagged in emerged_direction_support;\n     treat as low-confidence until reviewed. Set exclude_unsupported_emerged<-TRUE to drop.\n")
    } else {
      cat(sprintf("  Emerged-peak support check: all %d Emerged calls confirmed by a called peak (either tier) or the depth+magnitude rescue.\n", n_emerged))
    }
  }

  emerged_df <- data.frame(
    chr=as.character(seqnames(gr_em)), start=start(gr_em)-1L, end=end(gr_em),
    name=mcols(gr_em)$name, lfc_nc14b=NA_real_, padj_nc14b=NA_real_,
    direction_nc14b="Not differential at nc14b",
    lfc_nc14late=mcols(gr_em)$log2FC, padj_nc14late=mcols(gr_em)$padj,
    matched_nc14late=TRUE, fate="Emerged",
    direction_nc14late=ifelse(!is.na(mcols(gr_em)$log2FC) & mcols(gr_em)$log2FC>0,
                              lbl_grpA, lbl_grpB),
    temporal_grpA=emrg_ta, temporal_grpB=emrg_tb,
    temporal_grpA_nearby=emrg_ta_nearby, temporal_grpB_nearby=emrg_tb_nearby,  # [v14]
    temporal_label=case_when(
      emrg_ta & emrg_tb                   ~ paste("Both",grpA,"&",grpB,"change (overlapping)"),
      emrg_ta                             ~ paste(grpA,"changes over time (overlapping)"),
      emrg_tb                             ~ paste(grpB,"changes over time (overlapping)"),
      emrg_ta_nearby & emrg_tb_nearby     ~ paste("Both",grpA,"&",grpB,"change (nearby)"),
      emrg_ta_nearby                      ~ paste(grpA,"changes over time (nearby)"),
      emrg_tb_nearby                      ~ paste(grpB,"changes over time (nearby)"),
      TRUE                                 ~ "Neither changes over time"),
    temporal_lfc_grpA=NA_real_, temporal_lfc_grpB=NA_real_,
    temporal_dir_grpA=NA_character_, temporal_dir_grpB=NA_character_,
    resolution_mechanism=NA_character_,
    emerged_prior_nc14b_lfc=emerged_prior_lfc,
    emerged_confidence=emerged_confidence,
    emerged_direction_support=emerged_direction_support,
    stringsAsFactors=FALSE)
  fate_df <- bind_rows(fate_df, emerged_df)

  # Explicit exclusion-list rows (construct/driver-locus artifacts) are
  # dropped UNCONDITIONALLY -- regardless of exclude_unsupported_emerged --
  # since these are known artifacts, not a statistical judgment call you'd
  # ever want visible even in exploratory mode.
  if (has_called_peaks) {
    n_before_excl <- nrow(fate_df)
    is_excluded_locus <- !is.na(fate_df$emerged_direction_support) &
                          grepl("^EXCLUDED", fate_df$emerged_direction_support)
    if (any(is_excluded_locus)) {
      fn_exclloc <- file.path(out_contrast, "emerged_excluded_construct_artifact.csv")
      dir.create(out_contrast, recursive=TRUE, showWarnings=FALSE)
      write.csv(fate_df[is_excluded_locus, ], fn_exclloc, row.names=FALSE)
      cat("  Saved construct-artifact exclusion record:", fn_exclloc, "\n")
    }
    fate_df <- fate_df[!is_excluded_locus, ]
    cat(sprintf("  Dropped %d Emerged loci overlapping the explicit exclusion list (unconditional, regardless of exclude_unsupported_emerged)\n",
               n_before_excl - nrow(fate_df)))
  }

  if (exclude_unsupported_emerged && has_called_peaks) {
    n_before <- nrow(fate_df)
    is_unsupported <- !is.na(fate_df$emerged_direction_support) &
                       grepl("^UNSUPPORTED", fate_df$emerged_direction_support)
    if (any(is_unsupported)) {
      fn_excl <- file.path(out_contrast, "emerged_excluded_unsupported.csv")
      dir.create(out_contrast, recursive=TRUE, showWarnings=FALSE)
      write.csv(fate_df[is_unsupported, ], fn_excl, row.names=FALSE)
      cat("  Saved excluded-locus record:", fn_excl, "\n")
    }
    fate_df <- fate_df[!is_unsupported, ]
    cat(sprintf("  Dropped %d unsupported Emerged loci from fate_df (exclude_unsupported_emerged=TRUE)\n",
               n_before - nrow(fate_df)))
  }
}
fate_df <- as.data.frame(ungroup(fate_df))

# ── HCR-validated locus allowlist ──────────────────────────────────────────────
# The structural inverse of emerged_exclude_loci above: loci with independent,
# orthogonal validation (HCR) that fall short of the standard DAR magnitude/
# significance bar are added back explicitly, rather than by lowering
# log2fc_threshold/pvalue_threshold globally (which would relax the bar for
# every one of the ~185,000 tested windows in this dataset, not just the
# validated ones). This keeps statistical rigor intact for the dataset as a
# whole while still representing loci you have independent evidence for.
#
# Every allowlisted row carries a validation_source column (added to EVERY
# row in fate_df, not just these -- NA for normal DAR-called rows, populated
# only here) so downstream gene lists/figures/reviewers can immediately see
# and evaluate which rows rest on this basis versus standard DAR calling.
#
# Add entries as you validate more loci. Required columns: chr/start/end,
# direction_nc14b (use lbl_grpA or lbl_grpB, matching whichever genotype the
# HCR data shows as open), fate (typically "Maintained" for a consistent-
# direction, HCR-confirmed difference -- see below), lfc_nc14b/lfc_nc14late
# (best available point estimates, e.g. from whichever method/mode gave the
# lowest p-value even if it didn't clear filter_sig -- these are for display/
# magnitude-shift plots, not further statistical testing), and
# validation_source (free text -- cite the HCR experiment/figure).
hcr_validated_loci <- data.frame(
  chr   = c("chr3R"),
  start = c(13291398),   # EXTENDED200 merged window (matches the peak that
  end   = c(13294806),   # showed real signal in coverage tracks/HCR)
  direction_nc14b = c(lbl_grpB),   # BOTCv-open, per HCR and the consistent-
                                    # direction limma/DESeq2 point estimates
  fate  = c("Maintained"),         # same-direction at both timepoints per
                                    # HCR and available (sub-threshold) ATAC data
  lfc_nc14b    = c(-0.357040164562197),   # EXTENDED200_limma_BOTv_vs_BOTCv, p=0.0176
  lfc_nc14late = c(-0.822003471788508),   # NARROW_limma_BOTv_vs_BOTCv_nc14late, p=0.441
  matched_nc14late = c(TRUE),
  validation_source = c("HCR-validated (WntD); ATAC signal below standard DAR magnitude threshold (|log2FC|>=0.5) at both timepoints, and DESeq2/limma disagree sharply at nc14b (p=0.018 vs p=0.566) -- included on the basis of independent HCR confirmation, not ATAC statistical significance. See conversation/lab notes for HCR figure reference."),
  stringsAsFactors = FALSE
)
# TODO: fill in a proper HCR figure/experiment reference in validation_source
# above before this goes into any manuscript-facing output.

if (nrow(hcr_validated_loci) > 0) {
  fate_df$validation_source <- NA_character_
  fate_df <- dplyr::bind_rows(fate_df, hcr_validated_loci)
  cat(sprintf("  [HCR-ALLOWLIST] Added %d independently-validated locus/loci to fate_df (below standard DAR threshold, included on HCR evidence -- see validation_source column).\n",
             nrow(hcr_validated_loci)))
}

# ── Direction-split fate labels and colors ────────────────────────────────────
# Used by A6 (alluvial), A7 (stacked bar), and A9 (scatter) -- defined once,
# unconditionally, so all three plots draw from the same color source instead
# of three separately-maintained palettes drifting out of sync.
#
# "Deepened"/"Maintained"/"Reversed" alone don't say which genotype the peak
# ends up more open in. Color by the ENDING genotype so a color always means
# "this peak is BOTv-open" (green) or "this peak is BOTCv-open" (pink) at
# nc14late, consistently everywhere. Deepened/Maintained never change sign,
# so ending direction == starting direction (direction_nc14b) for those two;
# Reversed's ending direction is the opposite of direction_nc14b.
# ── PI-requested label scheme (literal strings, not derived from grpA/grpB,
# since the requested format "BOT.v"/"BOTC.v" doesn't match those variables'
# text -- if grpA/grpB genotype names ever change, these must be updated by
# hand). IMPORTANT: Reversed labels are keyed by ENDING genotype (same
# convention as the colors below, which were already correct) -- lbl_rev_A
# started grpA(BOTv)-open and ends grpB(BOTCv)-open, so it becomes
# "Reversed to BOTC.v", not "Reversed to BOT.v". Getting this backwards would
# recreate exactly the kind of direction-mislabeling bug fixed earlier for
# the Emerged category, just in a new place -- double-check against the
# "started grpA-open, ends grpB-open" comments before changing further.
lbl_emrg_A <- "Emerged BOT.v bias"
lbl_emrg_B <- "Emerged BOTC.v bias"
lbl_mnt_A  <- "Maintained BOT.v bias"     # ends grpA-open (also covers what
                                          # used to be split out as "Deepened")
lbl_mnt_B  <- "Maintained BOTC.v bias"    # ends grpB-open (same)
lbl_rev_A  <- "Reversed to BOTC.v"        # started grpA-open, ENDS grpB-open
lbl_rev_B  <- "Reversed to BOT.v"         # started grpB-open, ENDS grpA-open

fate_df$fate_split <- dplyr::case_when(
  fate_df$fate=="Maintained" & fate_df$direction_nc14b==lbl_grpA ~ lbl_mnt_A,
  fate_df$fate=="Maintained" & fate_df$direction_nc14b==lbl_grpB ~ lbl_mnt_B,
  fate_df$fate=="Reversed"   & fate_df$direction_nc14b==lbl_grpA ~ lbl_rev_A,
  fate_df$fate=="Reversed"   & fate_df$direction_nc14b==lbl_grpB ~ lbl_rev_B,
  TRUE ~ fate_df$fate)   # Converged / Emerged / Unknown pass through unchanged

# Master direction-split color/order map. Maintained-BOTCv keeps its former
# "Deepened" shade now that the two categories are merged (more visually
# distinct than the old lighter Maintained shade); everything BOTv-
# associated is a distinct green instead of inheriting the same pink.
fate_colors_bar <- c(
  setNames("#2ca02c", lbl_mnt_A),  # BOTv Maintained  -- bright green (former Deepened shade)
  setNames("#8B0057", lbl_rev_A),  # BOTv Reversed    -- dark magenta (ends BOTCv-open)
  setNames("#c9579a", lbl_mnt_B),  # BOTCv Maintained -- medium-deep pink (former Deepened shade)
  setNames("#1a7a3d", lbl_rev_B),  # BOTCv Reversed   -- dark green (ends BOTv-open)
  "Converged" = "#CCCCCC",
  "Emerged"  = "#F39C12"
)
fate_order_bar <- c(lbl_mnt_A, lbl_rev_A,
                    lbl_mnt_B, lbl_rev_B,
                    "Converged", "Emerged")

# ── Converged breakdown (moved out of the has_alluvial block so A13's
# ChIPseeker/gene-list grouping can use it even when ggalluvial isn't
# installed) ──────────────────────────────────────────────────────────────
# "Converged" collapses every locus whose BOTv-vs-BOTCv difference
# disappeared by nc14late into one bucket -- this splits it into 6
# sub-categories using temporal_dir_grpA/grpB (significance-filtered
# BOTv_temporal/BOTCv_temporal DAR calls, computed in A4b above).
fate_df$resolved_breakdown <- NA_character_
resolved_idx_plot <- which(fate_df$fate == "Converged" &
                           fate_df$direction_nc14b %in% c(lbl_grpA, lbl_grpB))
if (length(resolved_idx_plot) > 0) {
  fate_df$resolved_breakdown[resolved_idx_plot] <- with(
    fate_df[resolved_idx_plot,], dplyr::case_when(
      !is.na(temporal_dir_grpA) & !is.na(temporal_dir_grpB) ~ "Converged: both shift",
      temporal_dir_grpA=="opens"  & is.na(temporal_dir_grpB) ~ "Converged: BOTv opens",
      temporal_dir_grpA=="closes" & is.na(temporal_dir_grpB) ~ "Converged: BOTv closes",
      is.na(temporal_dir_grpA)    & temporal_dir_grpB=="opens"  ~ "Converged: BOTCv opens",
      is.na(temporal_dir_grpA)    & temporal_dir_grpB=="closes" ~ "Converged: BOTCv closes",
      TRUE ~ "Converged: neither genotype shows a temporal shift"))
}

# ── Unfiltered-trend fallback for "Converged: neither genotype shows a temporal shift" ─────────
# temporal_dir_grpA/grpB only reflect SIGNIFICANT BOTv_temporal/BOTCv_temporal
# DAR calls -- a locus with a real but modest/noisy trend that never cleared
# the significance threshold lands in "neither genotype shows a temporal shift" indistinguishable
# from one that's genuinely flat. If <stem>_FULL_results.txt (unfiltered
# DESeq2/limma results, every tested peak regardless of significance) exists
# for BOTv_temporal and/or BOTCv_temporal, use the raw log2FC direction to
# reclassify those rows -- this can only ever MOVE a row OUT of "no temporal
# signal" into a directional bucket, never the reverse, and never touches
# rows that already have a significant call.
#
# fate_df$mechanism_confidence records whether each Converged row's mechanism
# came from a significant DAR call or just a trending-but-not-significant
# unfiltered value, so this distinction is never silently lost -- it's
# exposed in peak_fate_data.csv and the per-group gene-list exports.
fate_df$mechanism_confidence <- NA_character_
fate_df$mechanism_confidence[resolved_idx_plot] <- "significant DAR"

trend_thresh <- 0.25   # much looser than significance-based calls -- distinguishes
                       # "some detectable direction" from noise-level log2FC
full_ta <- read_full_results("BOTv_temporal")
full_tb <- read_full_results("BOTCv_temporal")

if (is.null(full_ta) && is.null(full_tb)) {
  none_idx0 <- which(fate_df$resolved_breakdown == "Converged: neither genotype shows a temporal shift")
  fate_df$mechanism_confidence[none_idx0] <- "no signal (unfiltered data unavailable)"
  cat("\n  [NOTE] BOTv_temporal_FULL_results.txt / BOTCv_temporal_FULL_results.txt not found --\n",
      "  'Converged: neither genotype shows a temporal shift' reflects significance-filtered DAR calls only.\n",
      "  See write_full_results() in the Step1 export snippet to enable this fallback.\n")
} else {
  none_idx <- which(fate_df$resolved_breakdown == "Converged: neither genotype shows a temporal shift")
  if (length(none_idx) > 0) {
    gr_none <- GRanges(seqnames=fate_df$chr[none_idx],
                       ranges=IRanges(fate_df$start[none_idx]+1L, fate_df$end[none_idx]))
    trend_a <- rep(NA_real_, length(none_idx))
    trend_b <- rep(NA_real_, length(none_idx))
    if (!is.null(full_ta)) {
      gr_full_ta <- GRanges(seqnames=full_ta$chr,
                            ranges=IRanges(full_ta$start+1L, full_ta$end),
                            log2FC=full_ta$log2FoldChange)
      j_a <- any_overlap_join_mean(gr_none, gr_full_ta)
      if (nrow(j_a)>0) trend_a[j_a$i_q] <- j_a$lfc
    }
    if (!is.null(full_tb)) {
      gr_full_tb <- GRanges(seqnames=full_tb$chr,
                            ranges=IRanges(full_tb$start+1L, full_tb$end),
                            log2FC=full_tb$log2FoldChange)
      j_b <- any_overlap_join_mean(gr_none, gr_full_tb)
      if (nrow(j_b)>0) trend_b[j_b$i_q] <- j_b$lfc
    }
    a_hit <- !is.na(trend_a) & abs(trend_a) > trend_thresh
    b_hit <- !is.na(trend_b) & abs(trend_b) > trend_thresh
    new_breakdown <- dplyr::case_when(
      a_hit & b_hit                    ~ "Converged: both shift",
      a_hit & trend_a > 0               ~ "Converged: BOTv opens",
      a_hit & trend_a < 0               ~ "Converged: BOTv closes",
      b_hit & trend_b > 0               ~ "Converged: BOTCv opens",
      b_hit & trend_b < 0               ~ "Converged: BOTCv closes",
      TRUE                               ~ "Converged: neither genotype shows a temporal shift")
    new_confidence <- ifelse(new_breakdown == "Converged: neither genotype shows a temporal shift",
                             "no signal (checked unfiltered)",
                             "trending, not significant")
    fate_df$resolved_breakdown[none_idx]     <- new_breakdown
    fate_df$mechanism_confidence[none_idx]   <- new_confidence
    n_reclassified <- sum(new_breakdown != "Converged: neither genotype shows a temporal shift")
    cat(sprintf("\n  Unfiltered-trend fallback: reclassified %d/%d 'neither genotype shows a temporal shift' peaks using raw (non-significant) log2FC trend (|log2FC|>%.2f)\n",
               n_reclassified, length(none_idx), trend_thresh))
  }
}

# ── A6. Alluvial ──────────────────────────────────────────────────────────────

n_nc14b_total <- sum(fate_df$direction_nc14b %in% c(lbl_grpA, lbl_grpB))

alluvial_df <- fate_df %>%
  filter(direction_nc14b != "Unknown") %>%
  dplyr::count(direction_nc14b, fate, name="n") %>%
  group_by(direction_nc14b) %>%
  mutate(total=sum(n), pct=100*n/total) %>% ungroup() %>%
  mutate(direction_nc14b=factor(direction_nc14b,
                                levels=c(lbl_grpA,lbl_grpB,"Not differential at nc14b")),
         fate=factor(fate, levels=fate_order))

if (has_alluvial) {
  library(ggalluvial)

  # Alluvial-specific extension of the shared fate_colors_bar/fate_order_bar
  # map: adds the Emerged direction-split (Emerged has no "starting" state,
  # so it isn't part of fate_split / fate_colors_bar -- it's only split here,
  # locally, for the alluvial left-axis "De novo" source nodes).
  fate_colors_split <- fate_colors_bar
  fate_colors_split[lbl_emrg_A] <- "#4aab4a"
  fate_colors_split[lbl_emrg_B] <- "#d44fa8"
  fate_colors_split <- fate_colors_split[names(fate_colors_split) != "Emerged"]

  fate_order_split <- c(lbl_mnt_A,lbl_mnt_B,lbl_rev_A,lbl_rev_B,
                        "Converged", lbl_emrg_A, lbl_emrg_B)
  fate_order_nores <- c(lbl_mnt_A,lbl_mnt_B,lbl_rev_A,lbl_rev_B,
                        lbl_emrg_A, lbl_emrg_B)

  dir_colors_v2 <- c(setNames(col_grpA, lbl_grpA),
                     setNames(col_grpB, lbl_grpB),
                     "Not differential at nc14b" = "#F39C12")

  # Emerged direction split from fate_df (already correctly signed)
  emrg_rows  <- fate_df[fate_df$direction_nc14b == "Not differential at nc14b", ]
  n_emrg_A   <- sum(emrg_rows$direction_nc14late == lbl_grpA, na.rm=TRUE)
  n_emrg_B   <- sum(emrg_rows$direction_nc14late == lbl_grpB, na.rm=TRUE)
  n_botA_start <- sum(fate_df$direction_nc14b == lbl_grpA)
  n_botB_start <- sum(fate_df$direction_nc14b == lbl_grpB)
  n_resolved_v <- sum(fate_df$fate == "Converged" &
                      fate_df$direction_nc14b %in% c(lbl_grpA, lbl_grpB))

  base_alluv <- fate_df %>%
    dplyr::filter(direction_nc14b %in% c(lbl_grpA, lbl_grpB)) %>%
    dplyr::count(direction_nc14b, fate_split, name="n") %>%
    dplyr::rename(fate_base = fate_split)

  emrg_alluv <- data.frame(
    direction_nc14b = "Not differential at nc14b",
    fate_base = c(lbl_emrg_A, lbl_emrg_B),
    n = c(n_emrg_A, n_emrg_B), stringsAsFactors=FALSE)

  make_alluv_v2 <- function(excl_resolved=FALSE) {
    df <- if (excl_resolved) dplyr::filter(base_alluv, fate_base != "Converged") else base_alluv
    df %>% bind_rows(emrg_alluv) %>%
      dplyr::group_by(direction_nc14b) %>%
      dplyr::mutate(total=sum(n), pct=round(100*n/total,1)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        direction_nc14b = factor(direction_nc14b,
                                 levels=c(lbl_grpA, lbl_grpB, "Not differential at nc14b")),
        fate_base = factor(fate_base,
                           levels=if(excl_resolved) fate_order_nores else fate_order_split))
  }

  alluv_v2_full  <- make_alluv_v2(FALSE)
  alluv_v2_nores <- make_alluv_v2(TRUE)

  n_interesting_v   <- sum(fate_df$direction_nc14b %in% c(lbl_grpA,lbl_grpB) & fate_df$fate!="Converged")
  pct_A_nores <- round(100*sum(fate_df$direction_nc14b==lbl_grpA & fate_df$fate!="Converged")/n_botA_start,1)
  pct_B_nores <- round(100*sum(fate_df$direction_nc14b==lbl_grpB & fate_df$fate!="Converged")/n_botB_start,1)

  alluv_theme_v2 <- theme_minimal(base_size=11) +
    theme(panel.grid=element_blank(),
          axis.text.x=element_text(size=10,face="bold",colour="grey20"),
          axis.text.y=element_blank(), axis.ticks=element_blank(),
          plot.title=element_text(size=13,face="bold"),
          plot.subtitle=element_text(size=8.5,colour="grey40",lineheight=1.35),
          plot.caption=element_text(size=7.5,colour="grey55"),
          legend.position="none")

  # Master stratum color lookup used by all three alluvial plots
  # Keys must NOT contain \n — after_stat(stratum) strips \n before lookup
  stratum_colors_alluv <- c(
    setNames(col_grpA, lbl_grpA),
    setNames(col_grpB, lbl_grpB),
    "Not differential at nc14b" = "#F39C12",
    fate_colors_split
  )

  # ── Wide → lode form helper ───────────────────────────────────────────────
  # Wide format (axis1=/axis2=) mis-routes ribbons when multiple direction
  # groups share a fate stratum. Explicit lode form with integer alluvium IDs
  # guarantees each ribbon connects the correct (direction, fate) pair.
  wide_to_lode_v <- function(wide_df, dir_levels, fate_levels,
                              grey_resolved=FALSE) {
    wide_df <- wide_df %>%
      dplyr::mutate(
        alluvium   = dplyr::row_number(),
        direction_nc14b = factor(as.character(direction_nc14b), levels=dir_levels),
        fate_base  = factor(as.character(fate_base), levels=fate_levels),
        flow_alpha = dplyr::case_when(
          grey_resolved & as.character(fate_base)=="Converged" ~ 0.18,
          TRUE                                                  ~ 0.65
        )
      )
    dplyr::bind_rows(
      wide_df %>% dplyr::transmute(
        x=1L, stratum=as.character(direction_nc14b),
        alluvium=alluvium, y=n, flow_alpha=flow_alpha),
      wide_df %>% dplyr::transmute(
        x=2L, stratum=as.character(fate_base),
        alluvium=alluvium, y=n, flow_alpha=flow_alpha)
    ) %>%
      dplyr::mutate(
        stratum = factor(stratum, levels=c(dir_levels, fate_levels)),
        x       = factor(x, levels=c(1L, 2L))
      )
  }

  dir_levels_alluv <- c(lbl_grpA, lbl_grpB, "Not differential at nc14b")

  # Rebuild base_alluv with corrected De novo label (space, not \n)
  base_alluv <- fate_df %>%
    dplyr::filter(direction_nc14b %in% c(lbl_grpA, lbl_grpB)) %>%
    dplyr::count(direction_nc14b, fate_split, name="n") %>%
    dplyr::rename(fate_base = fate_split)

  emrg_alluv <- data.frame(
    direction_nc14b = "Not differential at nc14b",   # space not \n
    fate_base = c(lbl_emrg_A, lbl_emrg_B),
    n = c(n_emrg_A, n_emrg_B), stringsAsFactors=FALSE)

  wide_full_v <- dplyr::bind_rows(base_alluv, emrg_alluv)
  wide_nores_v <- dplyr::bind_rows(
    dplyr::filter(base_alluv, fate_base != "Converged"),
    emrg_alluv)

  lode_full_v  <- wide_to_lode_v(wide_full_v,  dir_levels_alluv, fate_order_split,
                                   grey_resolved=TRUE)
  lode_nores_v <- wide_to_lode_v(wide_nores_v, dir_levels_alluv, fate_order_nores,
                                   grey_resolved=FALSE)

  # ── PLOT 1: Full alluvial ─────────────────────────────────────────────────
  p_alluv_full <- ggplot(lode_full_v,
                         aes(x=x, stratum=stratum, alluvium=alluvium,
                             y=y, fill=stratum)) +
    # aes.flow="backward" → flows inherit LEFT (source) node color
    geom_flow(aes(alpha=flow_alpha), aes.flow="backward",
              width=0.38, knot.anchor=0.45, curve_type="sigmoid") +
    scale_alpha_identity() +
    geom_stratum(width=0.38, colour="white", linewidth=0.5) +
    scale_fill_manual(values=stratum_colors_alluv, guide="none") +
    geom_text(stat="stratum",
              aes(label=after_stat(paste0(stratum,"\n",prettyNum(count,big.mark=",")))),
              size=2.6, lineheight=1.25, colour="grey10") +
    scale_x_discrete(limits=c("1","2"),
                     labels=c("Direction at nc14b","Fate at nc14late"),
                     expand=c(0.20,0.20)) +
    labs(title=sprintf("Peak fate (Converged simplified): %s vs %s   nc14b \u2192 nc14late", grpA, grpB),
         subtitle=sprintf(
           "Starting: %s %s-open + %s %s-open  |  %s emerged de novo (%s opens in %s, %s in %s)\nConverged flows greyed \u2014 %s peaks (%.0f%% of starting pop) no longer a DAR at nc14late",
           prettyNum(n_botA_start,big.mark=","), grpA,
           prettyNum(n_botB_start,big.mark=","), grpB,
           prettyNum(n_emerged,big.mark=","),
           prettyNum(n_emrg_A,big.mark=","), grpA,
           prettyNum(n_emrg_B,big.mark=","), grpB,
           prettyNum(n_resolved_v,big.mark=","),
           100*n_resolved_v/(n_botA_start+n_botB_start)),
         y="Number of peaks", x=NULL,
         caption=sprintf("Maintained=same direction as nc14b  |  Reversed=direction flipped  |  Converged=no longer a DAR\nReversed colored by genotype it reverses INTO: %s Reversed (dark magenta, \u2192%s-open) vs %s Reversed (dark green, \u2192%s-open)  |  Emerged split by direction at nc14late: %s (light green) vs %s (light pink)",
                         grpA, grpB, grpB, grpA, grpA, grpB)) +
    alluv_theme_v2

  # ── PLOT 2: Active-fates zoom ──────────────────────────────────────────────
  p_alluv_nores <- ggplot(lode_nores_v,
                          aes(x=x, stratum=stratum, alluvium=alluvium,
                              y=y, fill=stratum)) +
    geom_flow(aes(alpha=flow_alpha), aes.flow="backward",
              width=0.38, knot.anchor=0.45, curve_type="sigmoid") +
    scale_alpha_identity() +
    geom_stratum(width=0.38, colour="white", linewidth=0.5) +
    scale_fill_manual(values=stratum_colors_alluv, guide="none") +
    geom_text(stat="stratum",
              aes(label=after_stat(paste0(stratum,"\nn=",prettyNum(count,big.mark=",")))),
              size=2.8, lineheight=1.25, colour="grey10") +
    scale_x_discrete(limits=c("1","2"),
                     labels=c("Direction at nc14b","Fate at nc14late"),
                     expand=c(0.22,0.22)) +
    labs(title=sprintf("Interesting fates only: %s vs %s   nc14b \u2192 nc14late", grpA, grpB),
         subtitle=sprintf(
           "Converged peaks excluded  |  Showing %s peaks (%s %s-open [%.1f%%], %s %s-open [%.1f%%])  +  %s emerged",
           prettyNum(n_interesting_v,big.mark=","),
           prettyNum(sum(fate_df$direction_nc14b==lbl_grpA & fate_df$fate!="Converged"),big.mark=","),
           grpA, pct_A_nores,
           prettyNum(sum(fate_df$direction_nc14b==lbl_grpB & fate_df$fate!="Converged"),big.mark=","),
           grpB, pct_B_nores,
           prettyNum(n_emerged,big.mark=",")),
         y="Number of peaks", x=NULL,
         caption=sprintf("Converged=%.0f%% of nc14b starting population \u2014 excluded to zoom in on active fate transitions\nReversed colored by genotype it reverses INTO: %s Reversed (dark magenta, \u2192%s-open) vs %s Reversed (dark green, \u2192%s-open)  |  Emerged split: opens in %s (light green) vs opens in %s (light pink)",
                         100*n_resolved_v/(n_botA_start+n_botB_start), grpA, grpB, grpB, grpA, grpA, grpB)) +
    alluv_theme_v2

  # ── PLOT 3: Panel ─────────────────────────────────────────────────────────
  p_alluv_panel <- (p_alluv_full  + labs(tag="A") + theme(plot.tag=element_text(size=14,face="bold"))) +
                   (p_alluv_nores + labs(tag="B") + theme(plot.tag=element_text(size=14,face="bold"))) +
    plot_annotation(
      title=sprintf("Peak fate: %s vs %s   nc14b \u2192 nc14late", grpA, grpB),
      subtitle="A: Full alluvial (Converged flows greyed)   |   B: Zoom \u2014 Converged excluded",
      theme=theme(plot.title=element_text(size=14,face="bold"),
                  plot.subtitle=element_text(size=9,colour="grey40")))

  # ── PLOT 4: Converged breakdown ────────────────────────────────────────────
  # "Converged" (93% of the starting population) collapses every locus whose
  # BOTv-vs-BOTCv difference disappeared by nc14late into one grey bucket --
  # but that disappearance happens for genuinely different reasons (BOTCv
  # catching up, BOTv fading, both moving, or no detectable temporal change
  # at all), and lumping them together hides which mechanism actually
  # dominates. fate_df$resolved_breakdown (computed unconditionally above,
  # before this has_alluvial block, so A13's ChIPseeker/gene-list grouping
  # can use it too) splits Converged into 6 sub-categories instead of 1.
  lbl_res <- c("Converged: BOTv opens","Converged: BOTv closes","Converged: both shift",
              "Converged: BOTCv opens","Converged: BOTCv closes","Converged: neither genotype shows a temporal shift")

  # Distinct, muted palette -- deliberately NOT reusing the bright green/pink
  # family from Deepened/Maintained/Reversed/Emerged, so "this used to be
  # Converged" stays visually distinguishable from "this is an active fate"
  # at a glance, even though some sub-categories lean green/pink thematically.
  fate_colors_resbreak <- c(
    fate_colors_split[setdiff(names(fate_colors_split),"Converged")],
    "Converged: BOTv opens"          = "#7fb3a0",  # soft teal-green
    "Converged: BOTv closes"         = "#a68a5b",  # muted olive/brown
    "Converged: both shift"          = "#8067b7",  # purple (distinct "combined" hue)
    "Converged: BOTCv opens"         = "#d99cc0",  # soft rose
    "Converged: BOTCv closes"        = "#8a6a7a",  # muted mauve
    "Converged: neither genotype shows a temporal shift"  = "#d9d9d9"   # light grey, distinguishable from stale "Converged" grey
  )
  fate_order_resbreak <- c(lbl_mnt_A,lbl_mnt_B,lbl_rev_A,lbl_rev_B,
                           lbl_res, lbl_emrg_A, lbl_emrg_B)
  stratum_colors_resbreak <- c(
    setNames(col_grpA, lbl_grpA),
    setNames(col_grpB, lbl_grpB),
    "Not differential at nc14b" = "#F39C12",
    fate_colors_resbreak
  )

  base_alluv_resbreak <- fate_df %>%
    dplyr::filter(direction_nc14b %in% c(lbl_grpA, lbl_grpB)) %>%
    dplyr::mutate(fate_final = ifelse(fate=="Converged", resolved_breakdown, fate_split)) %>%
    dplyr::count(direction_nc14b, fate_final, name="n") %>%
    dplyr::rename(fate_base = fate_final)

  wide_resbreak <- dplyr::bind_rows(base_alluv_resbreak, emrg_alluv)
  lode_resbreak <- wide_to_lode_v(wide_resbreak, dir_levels_alluv, fate_order_resbreak,
                                  grey_resolved=FALSE)

  n_res_botv_opens  <- sum(fate_df$resolved_breakdown=="Converged: BOTv opens", na.rm=TRUE)
  n_res_botv_closes <- sum(fate_df$resolved_breakdown=="Converged: BOTv closes", na.rm=TRUE)
  n_res_both        <- sum(fate_df$resolved_breakdown=="Converged: both shift", na.rm=TRUE)
  n_res_botcv_opens  <- sum(fate_df$resolved_breakdown=="Converged: BOTCv opens", na.rm=TRUE)
  n_res_botcv_closes <- sum(fate_df$resolved_breakdown=="Converged: BOTCv closes", na.rm=TRUE)
  n_res_none        <- sum(fate_df$resolved_breakdown=="Converged: neither genotype shows a temporal shift", na.rm=TRUE)

  p_alluv_resbreak <- ggplot(lode_resbreak,
                             aes(x=x, stratum=stratum, alluvium=alluvium,
                                 y=y, fill=stratum)) +
    geom_flow(aes(alpha=flow_alpha), aes.flow="backward",
              width=0.38, knot.anchor=0.45, curve_type="sigmoid") +
    scale_alpha_identity() +
    geom_stratum(width=0.38, colour="white", linewidth=0.5) +
    scale_fill_manual(values=stratum_colors_resbreak, guide="none") +
    geom_text(stat="stratum",
              aes(label=after_stat(paste0(stratum,"\n",prettyNum(count,big.mark=",")))),
              size=2.3, lineheight=1.2, colour="grey10") +
    scale_x_discrete(limits=c("1","2"),
                     labels=c("Direction at nc14b","Fate at nc14late (Converged broken down)"),
                     expand=c(0.20,0.20)) +
    labs(title=sprintf("Peak fate (Converged detailed \u2014 subcategories shown): %s vs %s   nc14b \u2192 nc14late", grpA, grpB),
         subtitle=sprintf(
           "Converged (%s peaks, %.0f%% of starting pop) split by which genotype's own accessibility moved during resolution\nBOTv opens=%s | BOTv closes=%s | both shift=%s | BOTCv opens=%s | BOTCv closes=%s | neither genotype shows a temporal shift=%s",
           prettyNum(n_resolved_v,big.mark=","),
           100*n_resolved_v/(n_botA_start+n_botB_start),
           prettyNum(n_res_botv_opens,big.mark=","), prettyNum(n_res_botv_closes,big.mark=","),
           prettyNum(n_res_both,big.mark=","),
           prettyNum(n_res_botcv_opens,big.mark=","), prettyNum(n_res_botcv_closes,big.mark=","),
           prettyNum(n_res_none,big.mark=",")),
         y="Number of peaks", x=NULL,
         caption="\"Converged\" alone only says the BOTv-vs-BOTCv difference disappeared -- these sub-categories say why.\n\"Both shift\" = both genotypes show a significant temporal DAR (any direction combination); \"neither genotype shows a temporal shift\" = neither genotype\nhas a detectable BOTv_temporal/BOTCv_temporal DAR at this locus (mechanism unresolvable without an unfiltered DESeq2 lookup).") +
    alluv_theme_v2 +
    theme(plot.subtitle=element_text(size=7.5,colour="grey40",lineheight=1.3))

  fn <- file.path(out_contrast,"peak_fate_alluvial_detailed_converged_subcategories.pdf")
  ggsave(fn, p_alluv_resbreak, width=12, height=8.5); cat("Saved:", fn, "\n")

  # ── PLOT 5: Converged breakdown, zoomed (no-temporal-signal excluded) ──────
  # Same idea as PLOT 2's zoom (drops Converged entirely) but less aggressive:
  # keeps every Converged sub-category that has SOME explanatory mechanism
  # (BOTv/BOTCv opens/closes, both shift) and only drops the one bucket that
  # genuinely has neither genotype shows a temporal shift in either genotype -- the "Converged" that
  # really is unexplained with current data, as opposed to the majority that
  # now have a concrete story.
  wide_resbreak_zoom <- dplyr::bind_rows(
    dplyr::filter(base_alluv_resbreak, fate_base != "Converged: neither genotype shows a temporal shift"),
    emrg_alluv)
  fate_order_resbreak_zoom <- setdiff(fate_order_resbreak, "Converged: neither genotype shows a temporal shift")
  lode_resbreak_zoom <- wide_to_lode_v(wide_resbreak_zoom, dir_levels_alluv,
                                       fate_order_resbreak_zoom, grey_resolved=FALSE)

  n_res_explained <- n_res_botv_opens + n_res_botv_closes + n_res_both +
                     n_res_botcv_opens + n_res_botcv_closes
  n_shown_zoom <- n_interesting_v + n_res_explained + n_emerged

  p_alluv_resbreak_zoom <- ggplot(lode_resbreak_zoom,
                                  aes(x=x, stratum=stratum, alluvium=alluvium,
                                      y=y, fill=stratum)) +
    geom_flow(aes(alpha=flow_alpha), aes.flow="backward",
              width=0.38, knot.anchor=0.45, curve_type="sigmoid") +
    scale_alpha_identity() +
    geom_stratum(width=0.38, colour="white", linewidth=0.5) +
    scale_fill_manual(values=stratum_colors_resbreak, guide="none") +
    geom_text(stat="stratum",
              aes(label=after_stat(paste0(stratum,"\nn=",prettyNum(count,big.mark=",")))),
              size=2.6, lineheight=1.2, colour="grey10") +
    scale_x_discrete(limits=c("1","2"),
                     labels=c("Direction at nc14b","Fate at nc14late (mechanism-explained only)"),
                     expand=c(0.20,0.20)) +
    labs(title=sprintf("Active + explained-Converged fates: %s vs %s   nc14b \u2192 nc14late", grpA, grpB),
         subtitle=sprintf(
           "\"Converged: neither genotype shows a temporal shift\" excluded (%s peaks, %.0f%% of starting pop) \u2014 showing %s peaks with an active fate or a concrete resolution mechanism",
           prettyNum(n_res_none,big.mark=","),
           100*n_res_none/(n_botA_start+n_botB_start),
           prettyNum(n_shown_zoom,big.mark=",")),
         y="Number of peaks", x=NULL,
         caption="Deepened/Maintained/Reversed = active BOTv-vs-BOTCv fate at nc14late (colored by ending genotype).\nConverged sub-categories = the BOTv-vs-BOTCv difference disappeared, but the genotype driving convergence is known.\nWhat's missing here (see peak_fate_alluvial_detailed_converged_subcategories.pdf for the full picture) is the majority of\nConverged peaks that have no detectable temporal DAR in either genotype.") +
    alluv_theme_v2 +
    theme(plot.subtitle=element_text(size=7.5,colour="grey40",lineheight=1.3))

  fn <- file.path(out_contrast,"peak_fate_alluvial_detailed_converged_subcategories_zoom.pdf")
  ggsave(fn, p_alluv_resbreak_zoom, width=12, height=8.5); cat("Saved:", fn, "\n")

  fn <- file.path(out_contrast,"peak_fate_alluvial_simplified_converged.pdf")
  ggsave(fn, p_alluv_full,  width=10, height=8); cat("Saved:", fn, "\n")
  fn <- file.path(out_contrast,"peak_fate_alluvial_active_fates.pdf")
  ggsave(fn, p_alluv_nores, width=9,  height=7); cat("Saved:", fn, "\n")
  fn <- file.path(out_contrast,"peak_fate_alluvial_panel.pdf")
  ggsave(fn, p_alluv_panel, width=19, height=8); cat("Saved:", fn, "\n")
}

# ── A7. Stacked bar ───────────────────────────────────────────────────────────

# Helper: wrap "More open in BOTv" → "More open in\nBOTv" (break before genotype)
wrap_dir_label <- function(x) sub(" (BOT)", "\n\\1", x, perl=TRUE)

# Direction-split version of alluvial_df (uses fate_split instead of fate) --
# needed for p_left/p_mid so Deepened/Maintained/Reversed each render in their
# ending-genotype color instead of one uniform color regardless of direction.
#
# tidyr::complete() guarantees every (direction_nc14b x fate_split) pair
# exists explicitly, filling n=0 for any combination that doesn't occur --
# without this, dplyr::count() only produces rows for combinations that
# actually have data, and ggplot's discrete scales drop unused factor levels
# by default. A rare category (e.g. n=1) could disappear from p_mid's y-axis
# entirely in a given run and reappear in another, which is confusing and
# looks like a rendering bug even when the underlying number is just small.
all_directions_split <- c(lbl_grpA, lbl_grpB, "Not differential at nc14b")
alluvial_df_split <- fate_df %>%
  filter(direction_nc14b != "Unknown") %>%
  dplyr::count(direction_nc14b, fate_split, name="n") %>%
  tidyr::complete(direction_nc14b=all_directions_split,
                  fate_split=fate_order_bar, fill=list(n=0)) %>%
  # Not differential at nc14b only ever has fate_split=="Emerged" (a peak that's
  # de novo has no Deepened/Maintained/Reversed/Converged fate to complete
  # against) -- drop the nonsensical zero-filled combinations tidyr::complete()
  # otherwise creates for that direction.
  filter(direction_nc14b != "Not differential at nc14b" | fate_split=="Emerged") %>%
  group_by(direction_nc14b) %>%
  mutate(total=sum(n), pct=100*n/total) %>% ungroup() %>%
  mutate(direction_nc14b=factor(direction_nc14b,
                                levels=c(lbl_grpA,lbl_grpB,"Not differential at nc14b")),
         fate_split=factor(fate_split, levels=fate_order_bar))

left_df <- alluvial_df_split %>%
  filter(direction_nc14b != "Not differential at nc14b") %>%
  mutate(fate_split=factor(fate_split,levels=rev(fate_order_bar)),
         direction_nc14b=factor(direction_nc14b,levels=c(lbl_grpA,lbl_grpB)))

# Categories too thin to label on the bar itself (<4% of their starting
# direction) -- called out in the caption instead, so a real (if small)
# number is never just silently invisible.
tiny_df <- left_df %>% filter(n>0, pct<4) %>%
  arrange(direction_nc14b, desc(pct)) %>%
  mutate(lbl=sprintf("%s: %s (n=%d, %.2f%%)",
                     gsub("\n"," ",wrap_dir_label(as.character(direction_nc14b))),
                     as.character(fate_split), n, pct))
tiny_caption <- if (nrow(tiny_df)>0)
  paste0("Too small to label on the bar (<4%): ", paste(tiny_df$lbl, collapse="  |  ")) else NULL

p_left <- ggplot(left_df, aes(x=direction_nc14b,y=pct,fill=fate_split)) +
  geom_bar(stat="identity",width=0.6,colour="white",linewidth=0.3) +
  geom_text(aes(label=ifelse(pct>=4,paste0(round(pct,1),"%"),"")),
            position=position_stack(vjust=0.5),size=2.8,colour="white",fontface="bold") +
  scale_fill_manual(values=fate_colors_bar,limits=fate_order_bar,name="Fate at nc14late") +
  scale_x_discrete(labels=wrap_dir_label) +
  scale_y_continuous(labels=percent_format(scale=1),expand=c(0,0),limits=c(0,102)) +
  labs(title="What happened\nto nc14b peaks?",x=NULL,y="% of starting direction",
       caption=tiny_caption) +
  theme_minimal(base_size=10) +
  theme(panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),
        plot.title=element_text(size=10,face="bold",lineheight=1.2),
        plot.caption=element_text(size=5.5,colour="grey40",hjust=0,lineheight=1.2),
        axis.text.x=element_text(size=9,colour=c(col_grpA,col_grpB),face="bold"),
        legend.position="right",legend.text=element_text(size=8),
        legend.title=element_text(size=8.5,face="bold"),legend.key.size=unit(0.8,"lines"))

mid_df <- alluvial_df_split %>% mutate(fate_split=factor(fate_split,levels=fate_order_bar))
p_mid <- ggplot(mid_df,aes(x=n,y=fate_split,colour=fate_split,shape=direction_nc14b)) +
  geom_segment(aes(x=0,xend=n,yend=fate_split),linewidth=0.6,colour="grey80") +
  geom_point(size=3.5,alpha=0.9) +
  geom_text(aes(label=prettyNum(n,big.mark=",")),hjust=-0.25,size=2.5,colour="grey30") +
  scale_colour_manual(values=fate_colors_bar,guide="none") +
  scale_y_discrete(drop=FALSE) +
  scale_shape_manual(values=setNames(c(16,17,15),c(lbl_grpA,lbl_grpB,"Not differential at nc14b")),
                     name="Starting direction",
                     labels=c(wrap_dir_label(lbl_grpA), wrap_dir_label(lbl_grpB),
                              "Not differential at nc14b")) +
  scale_x_continuous(expand=expansion(mult=c(0,0.25))) +
  labs(title="Absolute counts",x="Number of peaks",y=NULL) +
  theme_minimal(base_size=10) +
  theme(panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),
        plot.title=element_text(size=10,face="bold"),legend.position="bottom",
        legend.text=element_text(size=7.5),legend.key.size=unit(0.8,"lines"))

# NOTE: p_right intentionally stays on the UNSPLIT alluvial_df/fate. Its whole
# purpose is showing the BOTv/BOTCv composition WITHIN each fate category via
# the fill aesthetic -- switching its x-axis to fate_split would make most
# bars 100% single-color (since fate_split already resolves direction for
# Deepened/Maintained/Reversed), which would defeat the point of this panel.
right_df <- alluvial_df %>% mutate(fate=factor(fate,levels=fate_order)) %>%
  group_by(fate) %>% mutate(pct_of_fate=100*n/sum(n)) %>% ungroup()
p_right <- ggplot(right_df,aes(x=fate,y=pct_of_fate,fill=direction_nc14b)) +
  geom_bar(stat="identity",width=0.6,colour="white",linewidth=0.3) +
  geom_text(aes(label=ifelse(pct_of_fate>=5,paste0(round(pct_of_fate),"%"),"")),
            position=position_stack(vjust=0.5),size=2.7,colour="white",fontface="bold") +
  scale_fill_manual(values=c(setNames(col_grpA,lbl_grpA),setNames(col_grpB,lbl_grpB),
                              "Not differential at nc14b"="#F39C12"),
                    name="Starting direction",
                    labels=c(wrap_dir_label(lbl_grpA), wrap_dir_label(lbl_grpB),
                             "Not differential at nc14b")) +
  scale_x_discrete(labels=function(x) gsub(" ","\n",x)) +
  scale_y_continuous(labels=percent_format(scale=1),expand=c(0,0),limits=c(0,102)) +
  labs(title="Composition\nof each fate",x=NULL,y="% of fate class") +
  theme_minimal(base_size=10) +
  theme(panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),
        plot.title=element_text(size=10,face="bold",lineheight=1.2),
        legend.position="right",legend.text=element_text(size=8),
        legend.key.size=unit(0.8,"lines"))

p_bar <- (p_left|p_mid|p_right) +
  plot_annotation(
    title=sprintf("Peak fate: %s vs %s  nc14b → nc14late",grpA,grpB),
    subtitle=sprintf("nc14b: %s peaks  |  Emerged: %s peaks",
                     prettyNum(n_nc14b_total,big.mark=","),
                     prettyNum(n_emerged,big.mark=",")),
    caption="Left: fate proportions. Middle: absolute counts. Right: fate composition.",
    theme=theme(plot.title=element_text(size=13,face="bold"),
                plot.subtitle=element_text(size=8.5,colour="grey40"),
                plot.caption=element_text(size=7.5,colour="grey50")))
fn <- file.path(out_contrast,"peak_fate_stackedbar.pdf")
ggsave(fn, p_bar, width=13, height=6)
cat("Saved:", fn, "\n")

# ── A8. Temporal annotation plot ──────────────────────────────────────────────

has_temporal <- !is.null(df_temp_grpA) || !is.null(df_temp_grpB)
if (has_temporal) {
  temp_lvls <- c(paste(grpA,"changes over time"), paste(grpB,"changes over time"),
                 paste("Both",grpA,"&",grpB,"change"), "Neither changes over time")
  temp_colors <- c(setNames(col_grpA,paste(grpA,"changes over time")),
                   setNames(col_grpB,paste(grpB,"changes over time")),
                   setNames("#8B4513",paste("Both",grpA,"&",grpB,"change")),
                   "Neither changes over time"="#E0E0E0")
  temp_df <- fate_df %>%
    filter(direction_nc14b %in% c(lbl_grpA,lbl_grpB)) %>%
    dplyr::count(fate,direction_nc14b,temporal_label,name="n") %>%
    tidyr::complete(fate=setdiff(fate_order,"Emerged"), direction_nc14b=c(lbl_grpA,lbl_grpB),
                    temporal_label=temp_lvls, fill=list(n=0)) %>%
    group_by(fate,direction_nc14b) %>% mutate(pct=100*n/sum(n)) %>% ungroup() %>%
    mutate(fate=factor(fate,levels=setdiff(fate_order,"Emerged")),
           temporal_label=factor(temporal_label,levels=temp_lvls),
           direction_nc14b=factor(direction_nc14b,levels=c(lbl_grpA,lbl_grpB)),
           # Shorten strip labels for facets: "More open in BOTv" → "BOTv-open"
           dir_strip=factor(
             ifelse(direction_nc14b==lbl_grpA, paste0(grpA,"-open"), paste0(grpB,"-open")),
             levels=c(paste0(grpA,"-open"), paste0(grpB,"-open"))))
  p_temp <- ggplot(temp_df,aes(x=fate,y=pct,fill=temporal_label)) +
    geom_bar(stat="identity",width=0.65,colour="white",linewidth=0.3) +
    geom_text(aes(label=ifelse(pct>=6,paste0(round(pct),"%"),"")),
              position=position_stack(vjust=0.5),size=2.7,colour="white",fontface="bold") +
    facet_wrap(~dir_strip,ncol=2) +
    scale_fill_manual(values=temp_colors,name="Temporal co-change") +
    scale_x_discrete(labels=function(x) gsub(" ","\n",x)) +
    scale_y_continuous(labels=percent_format(scale=1),expand=c(0,0),limits=c(0,102)) +
    labs(title=sprintf("Temporal co-change by genotype fate (%s vs %s)",grpA,grpB),
         subtitle=paste0("Does the locus also change within ",grpA," or ",grpB," over time?"),
         x="Fate at nc14late",y="% of peaks in fate class") +
    theme_minimal(base_size=11) +
    theme(panel.grid.major.x=element_blank(),strip.text=element_text(size=10,face="bold"),
          plot.title=element_text(size=12,face="bold"),
          plot.subtitle=element_text(size=8.5,colour="grey40"),legend.position="right",
          legend.text=element_text(size=8.5),legend.key.size=unit(0.9,"lines"))
  fn <- file.path(out_contrast,"peak_fate_temporal_annot.pdf")
  ggsave(fn, p_temp, width=11, height=6)
  cat("Saved:", fn, "\n")
}

# ── A8b. Match confidence breakdown ────────────────────────────────────────────
# How much of each fate category comes from a directly significant nc14late
# DAR match (tier 1) vs. the trending-unfiltered (tier 2) or imputed-from-
# temporal-trajectories (tier 3) fallbacks built to recover peaks a single
# significance-filtered join would otherwise miss into "Converged". Lets you
# see at a glance which categories are mostly high-confidence direct
# measurements vs. mostly rescued/inferred, and decide how much to trust each
# one for downstream interpretation (e.g. gene lists, GO enrichment).
match_group_levels <- c(lbl_mnt_A, lbl_mnt_B, lbl_rev_A, lbl_rev_B,
                        paste0(grpA," Converged"), paste0(grpB," Converged"))
match_conf_levels <- c("significant DAR", "trending, not FDR-significant",
                       "imputed from BOTv/BOTCv temporal trajectories", "no nc14late match")
match_conf_colors <- c("significant DAR"="#2c7fb8",
                       "trending, not FDR-significant"="#7fcdbb",
                       "imputed from BOTv/BOTCv temporal trajectories"="#edf8b1",
                       "no nc14late match"="#d9d9d9")

mc_df <- fate_df %>%
  filter(direction_nc14b %in% c(lbl_grpA, lbl_grpB)) %>%
  mutate(match_group = ifelse(fate=="Converged",
                              paste0(gsub("More open in ","",direction_nc14b), " Converged"),
                              fate_split)) %>%
  dplyr::count(match_group, match_confidence, name="n") %>%
  tidyr::complete(match_group=match_group_levels, match_confidence=match_conf_levels,
                  fill=list(n=0)) %>%
  group_by(match_group) %>% mutate(total=sum(n), pct=100*n/total) %>% ungroup() %>%
  filter(total>0) %>%   # drop any category with genuinely zero peaks rather than show a blank/NaN bar
  mutate(match_group=factor(match_group,levels=rev(match_group_levels)),
         match_confidence=factor(match_confidence,levels=rev(match_conf_levels)))

if (nrow(mc_df)>0) {
  p_matchconf <- ggplot(mc_df, aes(x=match_group, y=pct, fill=match_confidence)) +
    geom_bar(stat="identity", width=0.7, colour="white", linewidth=0.3) +
    geom_text(aes(label=ifelse(pct>=6,paste0(round(pct),"%"),"")),
              position=position_stack(vjust=0.5), size=2.6, colour="grey15") +
    scale_fill_manual(values=match_conf_colors, limits=match_conf_levels, name="Match confidence") +
    scale_y_continuous(labels=percent_format(scale=1), expand=c(0,0), limits=c(0,102)) +
    coord_flip() +
    labs(title=sprintf("How was each fate call made? %s vs %s", grpA, grpB),
         subtitle="Significant nc14late DAR match vs. trending/imputed fallback tiers",
         x=NULL, y="% of peaks in category",
         caption="\"Converged\" here means matched_nc14late==FALSE (\"no nc14late match\") OR matched via a fallback tier\nthat still landed it in Converged -- see peak_fate_alluvial_detailed_converged_subcategories.pdf for the mechanism split.") +
    theme_minimal(base_size=11) +
    theme(panel.grid.major.y=element_blank(), legend.position="bottom",
          plot.title=element_text(size=12,face="bold"),
          plot.subtitle=element_text(size=8.5,colour="grey40"),
          plot.caption=element_text(size=7,colour="grey50"))

  fn <- file.path(out_contrast,"peak_fate_match_confidence.pdf")
  ggsave(fn, p_matchconf, width=9, height=6.5)
  cat("Saved:", fn, "\n")
}

# ── A9. LFC scatter ───────────────────────────────────────────────────────────

matched_df <- fate_df %>%
  filter(matched_nc14late, !is.na(lfc_nc14b), !is.na(lfc_nc14late),
         direction_nc14b %in% c(lbl_grpA,lbl_grpB)) %>%
  mutate(fate_split=factor(fate_split,levels=fate_order_bar))
if (nrow(matched_df)>0) {
  lim <- max(abs(c(matched_df$lfc_nc14b,matched_df$lfc_nc14late)),na.rm=TRUE)*1.05
  fate_in_scatter <- intersect(c(lbl_mnt_A,lbl_mnt_B,lbl_rev_A,lbl_rev_B),
                                levels(matched_df$fate_split))
  p_sc <- ggplot(matched_df,aes(x=lfc_nc14b,y=lfc_nc14late,colour=fate_split)) +
    geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey50",linewidth=0.7) +
    geom_hline(yintercept=0,colour="grey70",linewidth=0.4) +
    geom_vline(xintercept=0,colour="grey70",linewidth=0.4) +
    geom_point(alpha=0.6,size=1.4) +
    scale_colour_manual(values=fate_colors_bar,
                        breaks=fate_in_scatter,
                        name="Fate") +
    coord_fixed(xlim=c(-lim,lim),ylim=c(-lim,lim)) +
    labs(title=sprintf("LFC trajectory: %s vs %s   nc14b \u2192 nc14late",grpA,grpB),
         subtitle=sprintf("n=%d matched peaks",nrow(matched_df)),
         x=sprintf("log2FC at nc14b (+= more open in %s)",grpA),
         y=sprintf("log2FC at nc14late (+= more open in %s)",grpA),
         caption="Dashed diagonal: LFC unchanged. Color = fate, split by ending genotype (green=BOTv-open, pink=BOTCv-open).") +
    theme_minimal(base_size=11) +
    theme(panel.grid=element_line(colour="grey92"),
          plot.title=element_text(size=12,face="bold"),
          plot.subtitle=element_text(size=8.5,colour="grey40"))
  fn <- file.path(out_contrast,"peak_fate_lfc_scatter.pdf")
  ggsave(fn,p_sc,width=7,height=6.5); cat("Saved:",fn,"\n")
}

# ── A10. Magnitude shift ──────────────────────────────────────────────────────
# NOTE: this used to split Maintained vs Deepened at deepened_pct_thresh; now
# that those are merged into a single "Maintained" category (see A3), this
# shows the continuous magnitude-shift distribution for all Maintained peaks
# without a categorical split or the old (now-fictional) threshold line.
# pct_increase is still computed in A3 and still meaningful as a continuous
# value -- only the categorical cutoff that used to gate fate on it is gone.

mag_df <- fate_df %>%
  filter(fate=="Maintained", !is.na(lfc_nc14b), !is.na(lfc_nc14late),
         direction_nc14b %in% c(lbl_grpA,lbl_grpB)) %>%
  mutate(pct_increase_pct=100*pct_increase,   # pct_increase already computed in A3, as a fraction
         direction_nc14b=factor(direction_nc14b,levels=c(lbl_grpA,lbl_grpB)))
if (nrow(mag_df)>0) {
  p_mag <- ggplot(mag_df,
                  aes(x=direction_nc14b,y=pct_increase_pct,fill=direction_nc14b)) +
    geom_hline(yintercept=0,linetype="dashed",colour="grey50",linewidth=0.6) +
    geom_violin(alpha=0.7,width=0.8,colour=NA) +
    geom_boxplot(width=0.18,outlier.size=0.5,colour="grey30",fill="white") +
    scale_fill_manual(values=setNames(c(col_grpA,col_grpB),c(lbl_grpA,lbl_grpB)),guide="none") +
    labs(title="Magnitude shift for Maintained peaks",
         subtitle="% change in |fold-change| of the genotype-difference from nc14b to nc14late; positive = stronger difference, negative = weaker (but still same-direction, hence still Maintained)",
         x=NULL,y="% change in genotype-difference magnitude") +
    theme_minimal(base_size=11) +
    theme(panel.grid.major.x=element_blank(),
          plot.title=element_text(size=12,face="bold"))
  fn <- file.path(out_contrast,"peak_fate_magnitude_shift.pdf")
  ggsave(fn,p_mag,width=7,height=5); cat("Saved:",fn,"\n")
}

# ── A11. Peak widths ──────────────────────────────────────────────────────────

width_df <- fate_df %>%
  filter(direction_nc14b %in% c(lbl_grpA,lbl_grpB),
         fate %in% c("Maintained","Reversed","Converged")) %>%
  mutate(peak_width=end-start,
         fate=factor(fate,levels=c("Maintained","Reversed","Converged")),
         direction_nc14b=factor(direction_nc14b,levels=c(lbl_grpA,lbl_grpB)))
if (nrow(width_df)>0) {
  med_w <- width_df %>% group_by(fate,direction_nc14b) %>%
    summarise(med=median(peak_width,na.rm=TRUE),.groups="drop")
  p_wid <- ggplot(width_df,aes(x=fate,y=peak_width,fill=direction_nc14b)) +
    geom_violin(alpha=0.65,position=position_dodge(0.8),width=0.75,colour=NA) +
    geom_boxplot(width=0.15,position=position_dodge(0.8),
                 outlier.size=0.4,colour="grey30",fill="white") +
    geom_text(data=med_w,aes(y=med,label=paste0(round(med)," bp"),group=direction_nc14b),
              position=position_dodge(0.8),vjust=-0.5,size=2.6,colour="grey30") +
    scale_y_log10(labels=comma) +
    scale_fill_manual(values=setNames(c(col_grpA,col_grpB),c(lbl_grpA,lbl_grpB)),
                      name="Starting direction",
                      labels=c(wrap_dir_label(lbl_grpA), wrap_dir_label(lbl_grpB))) +
    scale_x_discrete(labels=function(x) gsub(" ","\n",x)) +
    labs(title="Peak width by fate and starting direction",
         subtitle="log10 y-axis | median annotated in bp",x=NULL,y="Peak width (bp, log10)") +
    theme_minimal(base_size=11) +
    theme(panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),
          plot.title=element_text(size=12,face="bold"))
  fn <- file.path(out_contrast,"peak_fate_widths.pdf")
  ggsave(fn,p_wid,width=8,height=5); cat("Saved:",fn,"\n")
}

# ── A12. Emerged direction ────────────────────────────────────────────────────

if (n_emerged>0) {
  emrg_df <- fate_df %>% filter(fate=="Emerged") %>%
    dplyr::count(direction_nc14late,name="n") %>%
    mutate(pct=100*n/sum(n),
           direction_nc14late=factor(direction_nc14late,levels=c(lbl_grpA,lbl_grpB)))
  p_em <- ggplot(emrg_df,aes(x=direction_nc14late,y=n,fill=direction_nc14late)) +
    geom_bar(stat="identity",width=0.55,colour="white") +
    geom_text(aes(label=sprintf("%s\n(%s%%)",prettyNum(n,big.mark=","),round(pct))),
              vjust=-0.3,size=3.5,fontface="bold") +
    scale_fill_manual(values=setNames(c(col_grpA,col_grpB),c(lbl_grpA,lbl_grpB)),guide="none") +
    scale_x_discrete(labels=function(x) gsub(" ","\n",x)) +
    scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
    labs(title=sprintf("Direction of %s de novo peaks at nc14late",
                       prettyNum(n_emerged,big.mark=",")),x=NULL,y="Number of peaks") +
    theme_minimal(base_size=11) +
    theme(panel.grid.major.x=element_blank(),plot.title=element_text(size=12,face="bold"))
  fn <- file.path(out_contrast,"peak_fate_emerged_direction.pdf")
  ggsave(fn,p_em,width=5,height=5); cat("Saved:",fn,"\n")
}

# ── A13. ChIPseeker annotation per fate ───────────────────────────────────────

if (has_chipseeker) {
  cat("\nRunning ChIPseeker annotation per fate group...\n")

  # Clean a resolved_breakdown value ("Converged: BOTv opens") into a group-
  # name-safe fragment ("BOTv_opens")
  clean_mechanism <- function(x) {
    x <- gsub("^Converged: ", "", x)
    gsub("[^A-Za-z0-9]+", "_", x)
  }

  # Build GRanges per group. Converged peaks (93% of the starting population)
  # get broken down into their resolved_breakdown mechanism sub-categories
  # instead of one undifferentiated "Converged" bucket -- e.g. "BOTv_Converged"
  # becomes "BOTv_Converged_BOTCv_opens", "BOTv_Converged_no_temporal_signal",
  # etc. Every other fate (Deepened/Maintained/Reversed/Emerged) is unaffected.
  #
  # Genotype token for the group/file name is derived explicitly here
  # ("BOTv"/"BOTCv"/"Emerged"), NOT by string-stripping direction_nc14b
  # (lbl_grpA/lbl_grpB). Those are the PI-facing DISPLAY labels and are
  # expected to keep changing wording over time (they already have, twice,
  # in this pipeline's history) -- deriving the internal/file-naming token
  # from them means every label wording change silently breaks every
  # downstream script that pattern-matches on gene-list filenames (e.g.
  # atac_rna_integration*.R's grepl("BOTv", fate_short) checks), exactly
  # what happened here when lbl_grpA changed from "More open in BOTv" to
  # "BOT.v-biased". Keeping this token hardcoded and independent of the
  # display label means display wording can change freely without breaking
  # any downstream consumer of these filenames.
  #
  # [v16] Emerged rows ALWAYS have direction_nc14b=="Not differential at
  # nc14b" (that's the definition of Emerged -- no nc14b DAR to have a
  # direction), so the case_when below could never match lbl_grpA/lbl_grpB
  # for them and every Emerged locus fell into the catch-all "Emerged" token
  # regardless of which genotype it actually opened in -- one undifferentiated
  # "Emerged_Emerged" group/file instead of a BOTv-open/BOTCv-open split.
  # The real per-locus direction for Emerged loci lives in direction_nc14late
  # (set explicitly in the emerged_df construction in A5), not
  # direction_nc14b -- added as an explicit fallback below, checked only for
  # fate=="Emerged" rows so it can't affect Converged/Maintained/Reversed
  # token assignment.
  fate_groups <- fate_df %>%
    filter(direction_nc14b %in% c(lbl_grpA, lbl_grpB, "Not differential at nc14b")) %>%
    mutate(geno_token = dplyr::case_when(
             direction_nc14b == lbl_grpA ~ "BOTv",
             direction_nc14b == lbl_grpB ~ "BOTCv",
             fate=="Emerged" & direction_nc14late==lbl_grpA ~ "BOTv",
             fate=="Emerged" & direction_nc14late==lbl_grpB ~ "BOTCv",
             TRUE                        ~ "Emerged"),   # true fallback: fate=="Emerged" but direction_nc14late somehow neither label (shouldn't happen; kept as a visible catch-all rather than silently mis-assigning)
           group = ifelse(
             fate=="Converged" & !is.na(resolved_breakdown),
             paste0(geno_token, "_Converged_", clean_mechanism(resolved_breakdown)),
             paste(geno_token, fate, sep="_")))

  ann_list_vent  <- list()
  gene_list_vent <- list()   # cleaned (deduped, SYMBOL-filtered) gene tables, for the Excel workbook
  for (grp in unique(fate_groups$group)) {
    sub_df <- fate_groups[fate_groups$group==grp, ]
    if (nrow(sub_df)<3) next
    gr_sub <- GRanges(seqnames=sub_df$chr,
                      ranges=IRanges(sub_df$start+1L, sub_df$end))
    cat(sprintf("  Annotating: %s (n=%d)...\n", grp, nrow(sub_df)))
    ann_df <- run_chipseeker(gr_sub, label=grp)
    if (!is.null(ann_df)) {
      ann_df$log2FC <- sub_df$lfc_nc14b[seq_len(nrow(ann_df))]
      ann_df$resolved_breakdown   <- sub_df$resolved_breakdown[seq_len(nrow(ann_df))]
      ann_df$mechanism_confidence <- sub_df$mechanism_confidence[seq_len(nrow(ann_df))]
      ann_df$match_confidence     <- sub_df$match_confidence[seq_len(nrow(ann_df))]
      ann_df$emerged_prior_nc14b_lfc <- sub_df$emerged_prior_nc14b_lfc[seq_len(nrow(ann_df))]
      ann_df$emerged_confidence      <- sub_df$emerged_confidence[seq_len(nrow(ann_df))]
      ann_list_vent[[grp]] <- ann_df
      gene_list_vent[[grp]] <- save_gene_list(ann_df, grp, out_gene_dir)
    }
  }

  if (length(ann_list_vent)>0) {
    p_ann <- plot_annotation_bars(ann_list_vent,
                                  title="Genomic feature annotation by fate group\n(BOTv vs BOTCv)")
    if (!is.null(p_ann)) {
      n_grps <- length(ann_list_vent)
      fn <- file.path(out_contrast,"peak_fate_chipseeker.pdf")
      ggsave(fn, p_ann, width=11, height=max(4, 1.2+n_grps*0.55))
      cat("Saved:", fn, "\n")
    }
  }

  # ── Combined multi-sheet Excel workbook: one sheet per group + a Summary ────
  if (has_excel && length(gene_list_vent)>0) {
    wb <- openxlsx::createWorkbook()

    n_peaks_per_group <- fate_groups %>% dplyr::count(group, name="n_peaks_total")
    summary_df <- data.frame(
      group   = names(gene_list_vent),
      n_genes = vapply(gene_list_vent, nrow, integer(1)),
      stringsAsFactors = FALSE
    ) %>% dplyr::left_join(n_peaks_per_group, by="group")
    openxlsx::addWorksheet(wb, "Summary")
    openxlsx::writeData(wb, "Summary", summary_df)
    openxlsx::setColWidths(wb, "Summary", cols=1:3, widths="auto")

    used_sheet_names <- character(0)
    for (grp in names(gene_list_vent)) {
      # Excel sheet names: <=31 chars, must be unique. Reserve space for the
      # dedup suffix BEFORE truncating the base name -- truncating to exactly
      # 31 first and THEN appending "_2" and re-truncating back to 31 just
      # cuts the suffix back off when the base is already at the limit,
      # leaving sheet_name unchanged and looping forever on long, similarly-
      # prefixed group names (exactly what the Converged-mechanism group
      # names produce, e.g. "BOTv_Converged_no_temporal_signal...").
      clean_name <- gsub("[^A-Za-z0-9_]","_", grp)
      sheet_name <- substr(clean_name, 1, 31)
      if (sheet_name %in% used_sheet_names) {
        suffix <- 2
        repeat {
          suffix_str <- paste0("_", suffix)
          candidate <- paste0(substr(clean_name, 1, 31 - nchar(suffix_str)), suffix_str)
          if (!(candidate %in% used_sheet_names)) { sheet_name <- candidate; break }
          suffix <- suffix + 1
        }
      }
      used_sheet_names <- c(used_sheet_names, sheet_name)

      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet_name, gene_list_vent[[grp]])
      openxlsx::setColWidths(wb, sheet_name,
                             cols=seq_len(ncol(gene_list_vent[[grp]])), widths="auto")
    }

    fn_xlsx <- file.path(out_gene_dir, "peak_fate_gene_lists_by_group.xlsx")
    openxlsx::saveWorkbook(wb, fn_xlsx, overwrite=TRUE)
    cat(sprintf("Saved: %s (%d sheets)\n", fn_xlsx, length(gene_list_vent)))
  } else if (!has_excel && length(gene_list_vent)>0) {
    cat("  [NOTE] openxlsx not installed -- skipping combined Excel workbook (install.packages('openxlsx'))\n")
  }
}

# ── A14. Statistical tests ────────────────────────────────────────────────────

cat("\nRunning statistical tests...\n")
stats_lines <- c("PART A: VENT STATISTICAL TESTS",
                 sprintf("Comparison: %s vs %s  nc14b → nc14late", grpA, grpB),
                 strrep("=",60), "")

# Fisher tests: Promoter and Distal Intergenic enrichment per fate
if (has_chipseeker && length(ann_list_vent)>0) {
  # Background = all annotated nc14b DARs
  all_gr <- GRanges(seqnames=fate_df$chr[fate_df$direction_nc14b %in% c(lbl_grpA,lbl_grpB)],
                    ranges=IRanges(fate_df$start[fate_df$direction_nc14b %in% c(lbl_grpA,lbl_grpB)]+1L,
                                   fate_df$end[fate_df$direction_nc14b %in% c(lbl_grpA,lbl_grpB)]))
  bg_ann <- run_chipseeker(all_gr, "background_vent")

  fisher_rows <- list()
  for (feat in c("Promoter","Distal Intergenic","Intron")) {
    for (grp in names(ann_list_vent)) {
      r <- fisher_feature_test(ann_list_vent[[grp]], feat, bg_ann, grp)
      if (!is.null(r)) fisher_rows[[paste(grp,feat)]] <- r
    }
  }
  if (length(fisher_rows)>0) {
    fisher_df <- do.call(rbind, fisher_rows)
    fisher_df$p_adj <- p.adjust(fisher_df$p_value, method="BH")
    fisher_df$sig   <- ifelse(fisher_df$p_adj<0.05,"*","")
    stats_lines <- c(stats_lines,
                     "Fisher's exact tests: feature enrichment vs background",
                     strrep("-",60),
                     capture.output(print(fisher_df, row.names=FALSE)), "")
    write.csv(fisher_df, file.path(out_contrast,"peak_fate_fisher_tests.csv"), row.names=FALSE)
    cat("  Saved Fisher test results\n")
  }
}

# Wilcoxon tests: LFC magnitude per fate
wilcox_rows <- list()
fate_df_dir <- fate_df %>% filter(direction_nc14b %in% c(lbl_grpA,lbl_grpB),
                                   !is.na(lfc_nc14b))
ref_fate <- "Converged"
for (f in setdiff(fate_order, ref_fate)) {
  sub  <- fate_df_dir$lfc_nc14b[fate_df_dir$fate==f]
  ref  <- fate_df_dir$lfc_nc14b[fate_df_dir$fate==ref_fate]
  if (length(sub)<3 || length(ref)<3) next
  wt <- wilcox.test(abs(sub), abs(ref), alternative="greater")
  wilcox_rows[[f]] <- data.frame(
    fate=f, vs=ref_fate,
    median_abs_lfc_fate=round(median(abs(sub)),3),
    median_abs_lfc_ref=round(median(abs(ref)),3),
    W=wt$statistic, p_value=signif(wt$p.value,3))
}
if (length(wilcox_rows)>0) {
  wilcox_df <- do.call(rbind,wilcox_rows)
  wilcox_df$p_adj <- p.adjust(wilcox_df$p_value, method="BH")
  wilcox_df$sig   <- ifelse(wilcox_df$p_adj<0.05,"*","")
  stats_lines <- c(stats_lines,
                   "Wilcoxon tests: |LFC| at nc14b compared to Converged baseline",
                   strrep("-",60),
                   capture.output(print(wilcox_df,row.names=FALSE)),"")
  write.csv(wilcox_df, file.path(out_contrast,"peak_fate_wilcoxon_tests.csv"), row.names=FALSE)
}

# LFC shift for maintained/deepened (nc14late vs nc14b)
if (nrow(mag_df)>0) {
  stats_lines <- c(stats_lines, "Magnitude shift (% change in |fold-change| of the genotype-difference, nc14b -> nc14late) by fate × direction")
  summ <- mag_df %>% group_by(fate,direction_nc14b) %>%
    summarise(n=n(), median_pct_increase=round(median(pct_increase_pct,na.rm=TRUE),1),
              p_vs_zero=signif(wilcox.test(pct_increase_pct)$p.value,3), .groups="drop")
  stats_lines <- c(stats_lines, capture.output(print(as.data.frame(summ),row.names=FALSE)),"")
}

fn_stats <- file.path(out_contrast,"peak_fate_stats_summary.txt")
writeLines(stats_lines, fn_stats)
cat("Saved:", fn_stats, "\n")

# ── A15. Data table ───────────────────────────────────────────────────────────

fn <- file.path(out_contrast,"peak_fate_data.csv")
write.csv(fate_df %>% arrange(direction_nc14b, fate,
                               desc(abs(lfc_nc14b)+abs(coalesce(lfc_nc14late,0)))),
          fn, row.names=FALSE)
cat("Saved:", fn, "\n")

################################################################################
# ── PART B: NON-VENT TEMPORAL ─────────────────────────────────────────────────
################################################################################

cat("\n", strrep("=",70), "\n", sep="")
cat("PART B: NON-VENT TEMPORAL PEAK FATE\n")
cat(strrep("=",70), "\n\n")

out_nv       <- file.path(out_dir, "peak_fate", "NonVent")
out_nv_genes <- file.path(out_nv, "gene_lists")
dir.create(out_nv_genes, recursive=TRUE, showWarnings=FALSE)

search_dirs_nv <- c(dar_dir, "./Temporal_Expansion_NonVent")

all_nv_rows <- do.call(rbind, lapply(nv_fate_pairs, function(fp) {
  df <- read_dar_full(fp$stem, search_dirs=search_dirs_nv)
  if (is.null(df)) return(NULL)
  df <- df[!is.na(df$log2FC), ]
  if (nrow(df)==0) return(NULL)

  anchor_df <- read_dar_full(fp$anchor_stem, search_dirs=search_dirs_nv)
  also_nc14b <- rep(FALSE, nrow(df))
  if (!is.null(anchor_df) && nrow(anchor_df)>0) {
    gr_t <- GRanges(seqnames=df$chr, ranges=IRanges(df$start+1L, df$end))
    gr_a <- GRanges(seqnames=anchor_df$chr, ranges=IRanges(anchor_df$start+1L, anchor_df$end))
    also_nc14b <- overlapsAny(gr_t, gr_a, minoverlap=1L)
  }

  df$also_nc14b <- also_nc14b
  df$direction  <- ifelse(df$log2FC>0, "Gains accessibility", "Loses accessibility")
  df$geno       <- fp$geno
  df$comparison <- sprintf("%s %s→%s", fp$geno, fp$tpA, fp$tpB)
  df$col_gain   <- fp$col_gain
  df$col_lose   <- fp$col_lose
  df
}))

if (!is.null(all_nv_rows) && nrow(all_nv_rows)>0) {

  # ── B1. Summary stacked bar ───────────────────────────────────────────────

  fate_summ_nv <- all_nv_rows %>%
    group_by(comparison, geno, direction, also_nc14b, col_gain, col_lose) %>%
    summarise(n=n(), .groups="drop") %>%
    group_by(comparison) %>%
    mutate(n_total=sum(n), pct=100*n/n_total) %>% ungroup() %>%
    mutate(fill_col=ifelse(direction=="Gains accessibility",col_gain,col_lose),
           nc14b_tag=ifelse(also_nc14b,"Also diff at nc14b","Temporal only"))

  p_nv <- ggplot(fate_summ_nv,
                 aes(x=comparison, y=pct, fill=fill_col,
                     alpha=ifelse(also_nc14b,0.65,1.0))) +
    geom_col(width=0.65,colour="white",linewidth=0.4,position="stack") +
    geom_hline(yintercept=50,linetype="dashed",colour="grey30",linewidth=0.7) +
    geom_text(aes(label=ifelse(pct>=6,sprintf("%.0f%%\nn=%d",pct,n),""),
                  colour=fill_col),
              position=position_stack(vjust=0.5),
              size=2.7,fontface="bold",lineheight=1.1) +
    scale_fill_identity() + scale_colour_identity() + scale_alpha_identity() +
    scale_y_continuous(labels=function(x) paste0(x,"%"),limits=c(0,105),expand=c(0,0)) +
    coord_flip() +
    labs(title="Non-vent temporal peak fate: gains vs losses",
         subtitle=paste0("Each bar = 100% of temporal DARs for that comparison\n",
                         "Solid = temporal-only  |  Semi-transparent = also genotype DAR at nc14b\n",
                         "Saturated = gaining accessibility over time; pale = losing"),
         x=NULL, y="% of temporal DARs") +
    theme_minimal(base_size=11) +
    theme(panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),
          axis.text.y=element_text(size=9.5,colour="grey15"),
          plot.title=element_text(size=12,face="bold"),
          plot.subtitle=element_text(size=8.5,colour="grey40",lineheight=1.35),
          plot.margin=margin(10,10,15,10))

  fn <- file.path(out_nv,"nonvent_temporal_fate.pdf")
  ggsave(fn, p_nv, width=10, height=max(4.5, 1.8+length(nv_fate_pairs)*0.9))
  cat("Saved:", fn, "\n")

  # ── B2. LFC distribution by comparison ──────────────────────────────────

  p_nv_lfc <- ggplot(all_nv_rows,
                     aes(x=log2FC, fill=direction, colour=direction)) +
    geom_density(alpha=0.45, linewidth=0.7) +
    geom_vline(xintercept=0, linetype="dashed", colour="grey40") +
    facet_wrap(~comparison, ncol=1, scales="free_y") +
    scale_fill_manual(values=c("Gains accessibility"="#145380",
                                "Loses accessibility"="#85C1E9"),guide="none") +
    scale_colour_manual(values=c("Gains accessibility"="#145380",
                                  "Loses accessibility"="#85C1E9"),guide="none") +
    labs(title="log2FC distribution for non-vent temporal DARs",
         subtitle="Positive = gains accessibility over time; negative = loses",
         x="log2FC (nc14late / nc14b)", y="Density") +
    theme_minimal(base_size=10) +
    theme(strip.text=element_text(size=9,face="bold"),
          panel.grid.minor=element_blank(),
          plot.title=element_text(size=11,face="bold"))
  fn <- file.path(out_nv,"nonvent_lfc_density.pdf")
  ggsave(fn, p_nv_lfc, width=7, height=3.5*length(nv_fate_pairs))
  cat("Saved:", fn, "\n")

  # ── B3. ChIPseeker + gene lists per non-vent group ───────────────────────

  if (has_chipseeker) {
    cat("\nRunning ChIPseeker annotation for non-vent groups...\n")
    ann_list_nv <- list()
    for (fp in nv_fate_pairs) {
      sub <- all_nv_rows[all_nv_rows$geno==fp$geno, ]
      if (nrow(sub)==0) next
      for (direc in unique(sub$direction)) {
        for (overlap in c(FALSE,TRUE)) {
          sub2 <- sub[sub$direction==direc & sub$also_nc14b==overlap, ]
          if (nrow(sub2)<3) next
          grp_nm <- sprintf("%s_%s_%s",
                            fp$geno,
                            ifelse(direc=="Gains accessibility","Gain","Lose"),
                            ifelse(overlap,"alsoNC14b","temporalOnly"))
          gr_sub <- GRanges(seqnames=sub2$chr, ranges=IRanges(sub2$start+1L,sub2$end))
          cat(sprintf("  Annotating: %s (n=%d)...\n", grp_nm, nrow(sub2)))
          ann_df <- run_chipseeker(gr_sub, label=grp_nm)
          if (!is.null(ann_df)) {
            ann_df$log2FC <- sub2$log2FC[seq_len(nrow(ann_df))]
            ann_list_nv[[grp_nm]] <- ann_df
            save_gene_list(ann_df, grp_nm, out_nv_genes)
          }
        }
      }
    }

    if (length(ann_list_nv)>0) {
      p_nv_ann <- plot_annotation_bars(ann_list_nv,
                                       title="Genomic annotation: non-vent temporal DARs")
      if (!is.null(p_nv_ann)) {
        fn <- file.path(out_nv,"nonvent_chipseeker.pdf")
        ggsave(fn, p_nv_ann, width=11, height=max(4, 1.2+length(ann_list_nv)*0.55))
        cat("Saved:", fn, "\n")
      }
    }
  }

  # ── B4. Statistical tests ────────────────────────────────────────────────

  cat("\nRunning non-vent statistical tests...\n")
  nv_stats <- c("PART B: NON-VENT TEMPORAL STATISTICAL TESTS", strrep("=",60), "")

  # Wilcoxon: |LFC| Gaining vs Losing per genotype
  for (fp in nv_fate_pairs) {
    sub <- all_nv_rows[all_nv_rows$geno==fp$geno, ]
    if (nrow(sub)<6) next
    gain <- abs(sub$log2FC[sub$direction=="Gains accessibility"])
    lose <- abs(sub$log2FC[sub$direction=="Loses accessibility"])
    if (length(gain)<3 || length(lose)<3) next
    wt <- wilcox.test(gain, lose, alternative="two.sided")
    nv_stats <- c(nv_stats,
      sprintf("%s: |LFC| Gain (median=%.2f, n=%d) vs Lose (median=%.2f, n=%d)  W=%.0f  p=%s",
              fp$geno, median(gain), length(gain), median(lose), length(lose),
              wt$statistic, signif(wt$p.value,3)))
  }
  nv_stats <- c(nv_stats, "")

  # Fisher: also-nc14b enrichment in Gain vs Lose per genotype
  nv_stats <- c(nv_stats, "Fisher: overlap with nc14b genotype DAR (Gain vs Lose)")
  for (fp in nv_fate_pairs) {
    sub <- all_nv_rows[all_nv_rows$geno==fp$geno, ]
    if (nrow(sub)<6) next
    n_gain_overlap <- sum(sub$direction=="Gains accessibility" & sub$also_nc14b)
    n_gain_only    <- sum(sub$direction=="Gains accessibility" & !sub$also_nc14b)
    n_lose_overlap <- sum(sub$direction=="Loses accessibility" & sub$also_nc14b)
    n_lose_only    <- sum(sub$direction=="Loses accessibility" & !sub$also_nc14b)
    mat <- matrix(c(n_gain_overlap, n_gain_only, n_lose_overlap, n_lose_only), 2, 2)
    ft <- tryCatch(fisher.test(mat), error=function(e) NULL)
    if (!is.null(ft))
      nv_stats <- c(nv_stats,
        sprintf("  %s: OR=%.2f  p=%s  (Gain:overlap=%d,only=%d | Lose:overlap=%d,only=%d)",
                fp$geno, ft$estimate, signif(ft$p.value,3),
                n_gain_overlap,n_gain_only,n_lose_overlap,n_lose_only))
  }
  nv_stats <- c(nv_stats,"")

  fn_nv_stats <- file.path(out_nv,"nonvent_stats_summary.txt")
  writeLines(nv_stats, fn_nv_stats)
  cat("Saved:", fn_nv_stats, "\n")

  # ── B5. Data table ────────────────────────────────────────────────────────

  fn <- file.path(out_nv,"nonvent_temporal_data.csv")
  write.csv(all_nv_rows %>%
              dplyr::select(chr,start,end,name,log2FC,geno,comparison,direction,also_nc14b) %>%
              arrange(geno,direction,desc(abs(log2FC))),
            fn, row.names=FALSE)
  cat("Saved:", fn, "\n")

} else {
  cat("  No non-vent temporal DAR files found — check stems and search dirs.\n")
}

################################################################################
# FINAL SUMMARY
################################################################################

cat("\n", strrep("=",70), "\n", sep="")
cat("COMPLETE. Output directories:\n")
cat("  Vent   :", out_contrast, "\n")
cat("  NonVent:", out_nv, "\n")
cat(strrep("=",70), "\n\n")
