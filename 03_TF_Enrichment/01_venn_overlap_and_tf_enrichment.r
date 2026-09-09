source("./00_tf_venn_shared_helpers.R")  # master palettes + recolored/label-safe Venn helpers + TF explanatory engine
################################################################################
# VENN SHARED REGIONS + TF ENRICHMENT — PER-CONDITION DAR SETS  (v3)
#
# Combined Venn + TF-enrichment analysis for one condition's DAR sets:
# specificity between the two focal genotypes, two background-correction
# methods (inferred vs empirical), and temporal (nc14b vs nc14late) dynamics.
# CONDITION_NAME/CONDITION_TAG below (currently: the ventralized BOTv/BOTCv
# split) are the only condition-specific inputs — point dar_base/temp_base at
# a different split's Step1 output and update dar_sets_s2 to reuse this same
# script structure for another condition.
#
# CANONICAL FILE — supersedes Venn_and_TF_Enrichment_Split2_Vent.r (v2).
#
# CHANGES IN v3:
#   - All Venn diagrams now recolored from the master genotype_colors /
#     concept_colors palette (00_tf_venn_shared_helpers.R) instead of ad hoc hex codes.
#     Two comparisons (nc14b vs nc14late; BOTCv above-inferred vs
#     above-empirical background) are conceptually about a *method* or
#     *timepoint* axis rather than a genotype axis, so they use the new
#     Timepoint_early/Timepoint_late and tolrm9_inferred/tolrm9_empirical/
#     Empirical_validation concept colors rather than genotype colors.
#   - FIX: category-name labels no longer get clipped at the plot edge —
#     same make_venn2()/make_venn3() label-wrapping + expanded coordinate
#     space + dynamic canvas width used in Split 1.
#   - TF panel switched from the old 9-TF placeholder list to the canonical
#     17-TF partner_peaks panel (PARTNER_PEAKS_PATHS), loaded via
#     load_partner_peaks() (rtracklayer::import() with read.table fallback).
#   - TF enrichment heatmap axis labels wrapped + margins expanded so DAR-set
#     and TF names no longer clip (same fix as Split 1).
#   - NEW Section B: TF explanatory breakdown (pies + bars + cross-set
#     summary) for every DAR set in dar_sets_s2, via the shared
#     run_explanatory_breakdown() / plot_explained_summary_across_sets()
#     engine — identical analysis to what Split 1 now produces.
################################################################################

CONDITION_TAG  <- "Split2_Vent"                    # short tag used in output paths
CONDITION_NAME <- "Split 2 Ventralized"            # human-readable label used in headers/titles

dar_base   <- "../Generate_fresh_counts/Output/SevenGeno_nc14b"
# FIX: temp_base was missing the "_v3" suffix -- it pointed at
# Split2_Vent_temporal (an old/incomplete directory), not
# Split2_Vent_temporal_v3 (the actual Step1 output location, confirmed by
# direct cp/wc -l verification in a prior debugging session). This directory
# holds BOTH the direct nc14b BOTv_vs_BOTCv comparison AND every temporal-
# derived stem (nc14late, BOTv_temporal, BOTCv_temporal, stable, emerging,
# resolving, temporal_divergent) -- so the missing suffix was silently
# shrinking every one of those sets, not just nc14late (was showing 1321,
# true value is 1371; nc14b BOTv_vs_BOTCv was showing 12,141 from dar_base,
# true value is 6,323, also sourced from here, not from dar_base).
temp_base  <- "../Generate_fresh_counts/Output/Split2_Vent_temporal_v3"
output_dir <- file.path(".", CONDITION_TAG, "Venn_Diagrams")
enrich_dir <- file.path(".", CONDITION_TAG, "TF_Enrichment")
explan_dir <- file.path(enrich_dir, "Explanatory_Breakdown")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(explan_dir, recursive = TRUE, showWarnings = FALSE)

peak_universe <- "../Enhanced_Run_wHets/ATAC_NarrowPeaks/SevenGeno_nc14b_union_peaks_FULL.bed"
# Confirmed real path (2026-06-21): byte-identical copy also exists at
# ../../Run_ChIPseq_12_04_25/ATAC_DAR_01.27.26/SevenGeno_nc14b_union_peaks_FULL.bed
# — same size/timestamp/line count (29,100 peaks), so either is fine; this
# one was chosen since it lives under the same parent project
# (United_peak_Universe_Code_10.28.25) as this script.

# ---------------------------------------------------------------------------
load_gr <- function(stem, base) {
  path <- file.path(base, paste0(stem, "_DARs.bed"))
  load_bed_generic(path)
}
overlap_count <- function(g1, g2) length(subsetByOverlaps(g1, g2, minoverlap = 1L))

################################################################################
# LOAD SETS
################################################################################

cat(sprintf("Loading DAR sets — %s...\n", CONDITION_NAME))
gr_AvsC_b   <- load_gr("BOTv_vs_BOTCv",                    temp_base)  # FIX: was dar_base (12,141, unfiltered SevenGeno_nc14b) -- BOTv_vs_BOTCv actually lives in Split2_Vent_temporal_v3 (true count 6,323), same dir as the temporal sets below
gr_AvsC_l   <- load_gr("BOTv_vs_BOTCv_nc14late",           temp_base)
gr_inf      <- load_gr("tolrm9_inferred_bg",                dar_base)
gr_emp      <- load_gr("tolrm9_empirical_bg",               dar_base)
gr_val      <- load_gr("tolrm9_validated_both_methods",     dar_base)
gr_Cv_inf   <- load_gr("BOTCv_above_inferred_bg",           dar_base)
gr_Cv_emp   <- load_gr("BOTCv_above_empirical_bg",          dar_base)
gr_Av_inf   <- load_gr("BOTv_above_inferred_bg",            dar_base)
gr_Cv_only  <- load_gr("BOTCv_only_crossvalidated",         dar_base)
gr_Av_only  <- load_gr("BOTv_only_crossvalidated",          dar_base)
gr_A_temp   <- load_gr("BOTv_temporal",                     temp_base)
gr_C_temp   <- load_gr("BOTCv_temporal",                    temp_base)
gr_stable   <- load_gr("BOTv_BOTCv_stable_both_timepoints", temp_base)
gr_emerging <- load_gr("BOTv_BOTCv_emerging_late_only",     temp_base)
gr_resolv   <- load_gr("BOTv_BOTCv_resolving_early_only",   temp_base)
gr_tdiv     <- load_gr("temporal_divergent_AvsC",           temp_base)
cat("\n")

################################################################################
# VENN DIAGRAMS — colors sourced from the master genotype_colors /
# concept_colors palette (see header for the per-Venn color rationale).
################################################################################

cat(sprintf("Generating Venn diagrams — %s...\n\n", CONDITION_NAME))

# 1. BOTv vs BOTCv at nc14b vs nc14late — this is a TIMEPOINT comparison of
#    the same BOTv-vs-BOTCv contrast, not a genotype comparison, so it uses
#    the Timepoint_early/Timepoint_late concept colors rather than the
#    BOTv/BOTCv genotype colors (which would misleadingly suggest the two
#    circles are different genotypes).
make_venn2(gr_AvsC_b, gr_AvsC_l,
           "BOTv/BOTCv nc14b", "BOTv/BOTCv nc14late",
           "BOTv_vs_BOTCv difference: nc14b vs nc14late",
           "Venn1_AvsC_early_vs_late.pdf",
           col1 = concept_colors[["Timepoint_early"]], col2 = concept_colors[["Timepoint_late"]])

# 2. tolrm9 background-correction method: inferred vs empirical
make_venn2(gr_inf, gr_emp,
           "tolrm9 inferred", "tolrm9 empirical",
           "tolrm9 background: inferred vs empirical",
           "Venn2_tolrm9_inferred_vs_empirical.pdf",
           col1 = concept_colors[["tolrm9_inferred"]], col2 = concept_colors[["tolrm9_empirical"]])

# 3. BOTCv specificity under the two background-correction methods
make_venn2(gr_Cv_inf, gr_Cv_emp,
           "BOTCv above inferred", "BOTCv above empirical",
           "BOTCv specificity: inferred vs empirical background subtraction",
           "Venn3_BOTCv_two_background_methods.pdf",
           col1 = genotype_colors[["BOTCv"]], col2 = concept_colors[["Empirical_validation"]])

# 4. Genotype specificity: BOTCv-only vs BOTv-only (cross-validated)
make_venn2(gr_Cv_only, gr_Av_only,
           "BOTCv-specific", "BOTv-specific",
           "Genotype specificity: BOTCv vs BOTv cross-validated",
           "Venn4_BOTCv_only_vs_BOTv_only.pdf",
           col1 = genotype_colors[["BOTCv"]], col2 = genotype_colors[["BOTv"]])

# 5. BOTv temporal change vs BOTCv temporal change — same genotype colors as
#    Venn4 for consistency, since this is still fundamentally a BOTv-vs-BOTCv
#    contrast, just measured as each genotype's own change over time.
make_venn2(gr_A_temp, gr_C_temp,
           "BOTv temporal", "BOTCv temporal",
           "Temporal changes: BOTv vs BOTCv",
           "Venn5_BOTv_vs_BOTCv_temporal.pdf",
           col1 = genotype_colors[["BOTv"]], col2 = genotype_colors[["BOTCv"]])

# 6. Three-way: stable / emerging / resolving temporal categories
if (length(gr_stable) > 0 && length(gr_emerging) > 0 && length(gr_resolv) > 0)
  make_venn3(gr_stable, gr_emerging, gr_resolv,
             "Stable", "Emerging", "Resolving",
             "Temporal categories: stable / emerging / resolving",
             "Venn6_temporal_categories_3way.pdf",
             col1 = concept_colors[["Stable"]], col2 = concept_colors[["Emerging"]],
             col3 = concept_colors[["Resolving"]])

# 7. Temporal-divergent DARs vs the static BOTv/BOTCv nc14b difference
make_venn2(gr_tdiv, gr_AvsC_b,
           "Temporal divergent", "BOTv/BOTCv nc14b",
           "Temporal divergent vs static BOTv/BOTCv difference",
           "Venn7_temporal_divergent_vs_static.pdf",
           col1 = concept_colors[["Temporal_divergent"]], col2 = genotype_colors[["BOTv"]])

################################################################################
# OVERLAP SUMMARY TABLE — with Jaccard and percent columns
################################################################################

cat("\nGenerating overlap summary table...\n")

ov_pairs <- list(
  list(g1 = gr_AvsC_b,  g2 = gr_AvsC_l,   l1 = "AvsC_nc14b",          l2 = "AvsC_nc14late"),
  list(g1 = gr_inf,     g2 = gr_emp,      l1 = "tolrm9_inferred",      l2 = "tolrm9_empirical"),
  list(g1 = gr_Cv_inf,  g2 = gr_Cv_emp,   l1 = "BOTCv_above_inferred", l2 = "BOTCv_above_empirical"),
  list(g1 = gr_Cv_only, g2 = gr_Av_only,  l1 = "BOTCv_only",           l2 = "BOTv_only"),
  list(g1 = gr_A_temp,  g2 = gr_C_temp,   l1 = "BOTv_temporal",        l2 = "BOTCv_temporal"),
  list(g1 = gr_stable,  g2 = gr_emerging, l1 = "Stable",                l2 = "Emerging"),
  list(g1 = gr_stable,  g2 = gr_resolv,   l1 = "Stable",                l2 = "Resolving"),
  list(g1 = gr_tdiv,    g2 = gr_AvsC_b,   l1 = "Temporal_divergent",    l2 = "AvsC_nc14b_static")
)

ov <- do.call(rbind, lapply(ov_pairs, function(p) {
  sh <- overlap_count(p$g1, p$g2)
  n1 <- length(p$g1); n2 <- length(p$g2)
  union_n <- n1 + n2 - sh
  data.frame(S1 = p$l1, N1 = n1, S2 = p$l2, N2 = n2, Shared = sh,
             Pct_of_S1 = round(100 * sh / max(n1, 1), 1),
             Pct_of_S2 = round(100 * sh / max(n2, 1), 1),
             Jaccard   = round(sh / max(union_n, 1), 3))
}))

write.table(ov, file.path(output_dir, "overlap_summary_Split2.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(ov)
cat("\nVenn diagrams saved to:", output_dir, "\n\n")

################################################################################
# TF ENRICHMENT — PER-CONDITION DAR SETS
#
# A. STATISTICAL ENRICHMENT — Fisher's exact test against the peak universe
# B. EXPLANATORY BREAKDOWN (NEW) — descriptive pies/bars: what % of DARs are
#    touched by partner TF binding, single-TF vs combinatorial, and which
#    TFs carry the most weight
################################################################################

dar_sets_s2 <- list(
  # --- "direct" genotype comparisons (BOTv vs BOTCv itself) + standalone
  # temporal sets — these already get Venn diagrams above (Venn1, Venn5)
  # but were missing from this list, so they never got a Fisher enrichment
  # row or an explanatory-breakdown panel. Added so both analyses cover
  # them too.
  list(stem = "BOTv_vs_BOTCv",                      dir = temp_base, label = "BOTv vs BOTCv (direct, nc14b)",    color = concept_colors[["Timepoint_early"]]),  # FIX: was dar_base (12,141 unfiltered) -- true set lives in temp_base (Split2_Vent_temporal_v3)
  list(stem = "BOTv_vs_BOTCv_nc14late",              dir = temp_base, label = "BOTv vs BOTCv (direct, nc14late)", color = concept_colors[["Timepoint_late"]]),
  list(stem = "BOTv_temporal",                       dir = temp_base, label = "BOTv temporal",                    color = genotype_colors[["BOTv"]]),
  list(stem = "BOTCv_temporal",                      dir = temp_base, label = "BOTCv temporal",                   color = genotype_colors[["BOTCv"]]),
  list(stem = "BOTCv_above_inferred_bg",           dir = dar_base,  label = "BOTCv above inferred bg",  color = genotype_colors[["BOTCv"]]),
  list(stem = "BOTCv_above_empirical_bg",          dir = dar_base,  label = "BOTCv above empirical bg", color = concept_colors[["Empirical_validation"]]),
  list(stem = "BOTCv_only_crossvalidated",         dir = dar_base,  label = "BOTCv-specific",            color = genotype_colors[["BOTCv_late"]]),
  list(stem = "BOTv_above_inferred_bg",            dir = dar_base,  label = "BOTv above inferred bg",    color = genotype_colors[["BOTv"]]),
  list(stem = "BOTv_only_crossvalidated",          dir = dar_base,  label = "BOTv-specific",              color = genotype_colors[["BOTv_late"]]),
  list(stem = "BOTv_BOTCv_stable_both_timepoints", dir = temp_base, label = "Stable (both TP)",           color = concept_colors[["Stable"]]),
  list(stem = "BOTv_BOTCv_emerging_late_only",     dir = temp_base, label = "Emerging (late only)",       color = concept_colors[["Emerging"]]),
  list(stem = "temporal_divergent_AvsC",           dir = temp_base, label = "Temporal divergent",         color = concept_colors[["Temporal_divergent"]])
)

# ---------------------------------------------------------------------------
fisher_enrich <- function(dar_gr, tf_gr, universe_gr) {
  if (length(dar_gr) == 0 || length(tf_gr) == 0 || length(universe_gr) == 0) return(NULL)
  in_dar  <- countOverlaps(dar_gr,      tf_gr, minoverlap = 1L) > 0
  in_univ <- countOverlaps(universe_gr, tf_gr, minoverlap = 1L) > 0
  a <- sum(in_dar)
  b <- sum(!in_dar)
  c <- sum(in_univ) - a
  d <- max(length(universe_gr) - a - b - c, 0)
  ft <- tryCatch(fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater"),
                 error = function(e) list(p.value = NA_real_, estimate = NA_real_))
  data.frame(
    n_dar                = length(dar_gr),
    n_tf_in_dar          = a,
    pct_dar_with_tf      = round(100 * a / max(length(dar_gr), 1), 1),
    pct_universe_with_tf = round(100 * sum(in_univ) / max(length(universe_gr), 1), 1),
    odds_ratio           = as.numeric(ft$estimate),
    pvalue               = ft$p.value
  )
}

cat("Loading peak universe...\n")
universe_gr <- tryCatch(load_bed_generic(peak_universe),
                        error = function(e) { message("Universe not found."); GRanges() })

cat("Loading partner TF peak files (canonical 17-TF panel)...\n")
tf_grs <- load_partner_peaks()
cat("\n")

################################################################################
# A. STATISTICAL ENRICHMENT
################################################################################

cat(sprintf("Running TF enrichment — %s...\n\n", CONDITION_NAME))

all_res <- list()
for (ds in dar_sets_s2) {
  dar_gr <- load_gr(ds$stem, ds$dir)
  if (length(dar_gr) == 0) next
  cat("  Testing:", ds$stem, "(", length(dar_gr), "DARs )\n")
  res <- do.call(rbind, lapply(names(tf_grs), function(tf_nm) {
    r <- fisher_enrich(dar_gr, tf_grs[[tf_nm]], universe_gr)
    if (is.null(r)) return(NULL)
    cbind(data.frame(dar_set = ds$stem, label = ds$label, tf = tf_nm), r)
  }))
  if (!is.null(res) && nrow(res) > 0) {
    res$padj <- p.adjust(res$pvalue, method = "BH")
    res$sig  <- factor(res$padj < 0.05, levels = c(FALSE, TRUE))
    all_res[[ds$stem]] <- res
    write.table(res, file.path(enrich_dir, paste0("TF_enrichment_", ds$stem, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  rm(dar_gr, res); gc(verbose = FALSE)
}

if (length(all_res) == 0) {
  # IMPORTANT: never call quit() from inside a sourced analysis script — it
  # kills the whole live R process, which RStudio reports as "R Session
  # Aborted" / "fatal error", even though no actual crash occurred. Skip
  # gracefully instead and print diagnostics so the real cause is visible.
  cat("\nNo enrichment results for any DAR set. Likely causes:\n")
  cat("  - tf_grs has", length(tf_grs), "TF(s) loaded",
      if (length(tf_grs) == 0) "  <-- 0 TFs loaded, check PARTNER_PEAKS_PATHS / working directory" else "", "\n")
  cat("  - universe_gr has", length(universe_gr), "peak(s) loaded",
      if (length(universe_gr) == 0) "  <-- 0 peaks loaded, check peak_universe path" else "", "\n")
  cat("  - current working directory:", getwd(), "\n")
  cat("  - peak_universe path tested:", peak_universe, "\n")
  cat("  - peak_universe file exists at that path?", file.exists(peak_universe), "\n")
  cat("\nSkipping heatmap section; continuing to Part B (explanatory breakdown).\n\n")
}

################################################################################
# SUMMARY HEATMAP — all sets x all TFs
# FIX: axis text wrapped + expanded margins so neither the DAR-set names
# (x-axis) nor the TF names (y-axis) get clipped, and canvas width scales
# with both the number of sets and the longest wrapped label.
################################################################################

if (length(all_res) > 0) {

  cat("\nGenerating enrichment heatmap...\n")

  combined <- do.call(rbind, all_res)
  combined$label <- factor(combined$label, levels = sapply(dar_sets_s2, `[[`, "label"))

  p_heat <- ggplot(combined, aes(x = label, y = tf,
                                  fill = log2(pmax(odds_ratio, 0.01)),
                                  size = sig)) +
    geom_point(shape = 21, colour = "white") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d6604d",
                         midpoint = 0, name = "log2(OR)") +
    scale_size_manual(values = c("FALSE" = 3, "TRUE" = 6), guide = "none") +
    scale_x_discrete(labels = function(x) wrap_label(x, 14)) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
          axis.text.y = element_text(size = 9),
          panel.grid  = element_blank(),
          plot.margin = margin(10, 20, 10, 10)) +
    labs(title    = sprintf("TF enrichment — %s", CONDITION_NAME),
         subtitle = "BOTCv and BOTv specificity layers + temporal categories\nCircle size = significant (BH padj<0.05) | Color = log2 odds ratio",
         x = NULL, y = "TF")

  ggsave(file.path(enrich_dir, "TF_enrichment_Split2_heatmap.pdf"), p_heat,
         width  = max(8, length(dar_sets_s2) * 1.1),
         height = max(5, length(tf_grs) * 0.45 + 2))
  cat("  Saved: TF_enrichment_Split2_heatmap.pdf\n")

  write.table(combined, file.path(enrich_dir, "TF_enrichment_Split2_all.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("TF enrichment heatmap saved.\n\n")
}

################################################################################
# B. TF EXPLANATORY BREAKDOWN (NEW) — what % of DARs in each set are touched
# by >=1 partner TF, and how is that explanatory power distributed across
# single-TF vs combinatorial binding and across individual TFs?
################################################################################

cat("Running TF explanatory breakdown (pies + bars)...\n\n")

explan_summaries <- list()
for (ds in dar_sets_s2) {
  dar_gr <- load_gr(ds$stem, ds$dir)
  if (length(dar_gr) == 0) next
  cat("  Explanatory breakdown:", ds$label, "\n")
  row <- run_explanatory_breakdown(dar_gr, tf_grs, label = ds$label,
                                    out_dir = explan_dir, file_stub = ds$stem)
  if (!is.null(row)) explan_summaries[[ds$stem]] <- row
  rm(dar_gr, row); gc(verbose = FALSE)
}

if (length(explan_summaries) > 0) {
  explan_df <- do.call(rbind, explan_summaries)
  write.table(explan_df, file.path(explan_dir, "Explanatory_summary_Split2.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  plot_explained_summary_across_sets(
    explan_df, explan_dir,
    fname     = "Explanatory_summary_across_sets_Split2.pdf",
    title_str = "TF explanatory power across ventralized DAR sets")
  cat("\n  Explanatory breakdown saved to:", explan_dir, "\n")
  print(explan_df)
} else {
  cat("  No explanatory breakdown results — check DAR sets and TF panel.\n")
}

cat(sprintf("\n%s Venn + TF enrichment + explanatory breakdown complete.\n\n", CONDITION_NAME))
