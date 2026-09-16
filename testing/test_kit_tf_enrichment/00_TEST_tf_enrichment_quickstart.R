################################################################################
# 00_TEST_tf_enrichment_quickstart.R
#
# Minimal reviewer-facing test for 00_tf_venn_shared_helpers.R.
#
# Runs the shared TF-enrichment / Venn helper functions on a small synthetic
# dataset (30 + 25 fake "DAR" peaks, 4 fake "TF" peak sets) so a reviewer can
# confirm the analysis functions work correctly without needing the real
# multi-gigabyte ATAC-seq/ChIP-seq data.
#
# WHAT THIS PROVES:
#   - load_bed_generic() correctly parses BED files into GRanges
#   - tf_explanatory_summary() correctly computes DAR x TF overlap stats
#   - make_venn2() produces a Venn diagram PDF with correct set sizes
#   - run_explanatory_breakdown() produces the full plot panel + tables
#
# EXPECTED OUTPUT (see test_data/README.md for how these numbers were derived):
#   Runt_specific DARs (n=30): 20 explained (66.7%), 15 single-TF, 5 multi-TF
#   Cic_specific  DARs (n=25): 22 explained (88.0%), 12 single-TF, 10 multi-TF
#   Runt_specific vs Cic_specific shared peaks (Venn): 10
#
# USAGE: Rscript 00_TEST_tf_enrichment_quickstart.R
# (run from inside this test_kit/ folder — paths below are relative to it)
################################################################################

# --- locate the real helper script -------------------------------------
# Adjust this path if you place the test_kit folder somewhere else relative
# to the repo. Defaults to sibling of ATACseq_Processing/.
helpers_path <- "../ATACseq_Processing/03_TF_Enrichment/00_tf_venn_shared_helpers.R"
if (!file.exists(helpers_path)) {
  stop("Can't find 00_tf_venn_shared_helpers.R at: ", helpers_path,
       "\nEdit `helpers_path` at the top of this script to point at it.")
}

output_dir <- "./test_output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source(helpers_path)

cat("\n================ TEST DATA SETUP ================\n")

# Override the real 17-TF panel with our 4 toy TFs for this test only.
TEST_PARTNER_PEAKS_PATHS <- list(
  TF_Runt = "./test_data/tf_partners/TF_Runt.bed",
  TF_Cic  = "./test_data/tf_partners/TF_Cic.bed",
  TF_Zld  = "./test_data/tf_partners/TF_Zld.bed",
  TF_Bcd  = "./test_data/tf_partners/TF_Bcd.bed"
)

cat("Loading toy TF partner peaks...\n")
tf_grs <- load_partner_peaks(TEST_PARTNER_PEAKS_PATHS)

cat("\nLoading toy DAR sets...\n")
dar_runt <- load_bed_generic("./test_data/DAR_Runt_specific.bed")
dar_cic  <- load_bed_generic("./test_data/DAR_Cic_specific.bed")
cat(sprintf("  Runt_specific DARs: %d peaks\n", length(dar_runt)))
cat(sprintf("  Cic_specific  DARs: %d peaks\n", length(dar_cic)))

cat("\n================ TEST 1: tf_explanatory_summary() ================\n")
summ_runt <- tf_explanatory_summary(dar_runt, tf_grs, label = "Runt_specific")
summ_cic  <- tf_explanatory_summary(dar_cic,  tf_grs, label = "Cic_specific")

cat(sprintf("Runt_specific: n_total=%d  explained=%d (%.1f%%)  single=%d  multi=%d\n",
            summ_runt$n_total, summ_runt$n_explained,
            100 * summ_runt$n_explained / summ_runt$n_total,
            summ_runt$n_single, summ_runt$n_multi))
cat(sprintf("Cic_specific:  n_total=%d  explained=%d (%.1f%%)  single=%d  multi=%d\n",
            summ_cic$n_total, summ_cic$n_explained,
            100 * summ_cic$n_explained / summ_cic$n_total,
            summ_cic$n_single, summ_cic$n_multi))

# --- sanity checks against hand-computed expected values ----------------
stopifnot(
  "Runt_specific n_total should be 30"    = summ_runt$n_total == 30,
  "Runt_specific n_explained should be 20" = summ_runt$n_explained == 20,
  "Cic_specific n_total should be 25"      = summ_cic$n_total == 25,
  "Cic_specific n_explained should be 22"  = summ_cic$n_explained == 22
)
cat("\n[PASS] tf_explanatory_summary() output matches expected values.\n")

cat("\n================ TEST 2: run_explanatory_breakdown() (full plots) ================\n")
breakdown_summary <- rbind(
  run_explanatory_breakdown(dar_runt, tf_grs, "Runt_specific", output_dir),
  run_explanatory_breakdown(dar_cic,  tf_grs, "Cic_specific",  output_dir)
)
print(breakdown_summary)
cat("[PASS] Explanatory breakdown plots written to ", output_dir, "\n")

cat("\n================ TEST 3: make_venn2() ================\n")
shared_n <- length(subsetByOverlaps(dar_runt, dar_cic, minoverlap = 1L))
cat(sprintf("Computed shared peaks (Runt_specific vs Cic_specific): %d (expected: 10)\n", shared_n))
stopifnot("Shared peak count should be 10" = shared_n == 10)

make_venn2(dar_runt, dar_cic, "Runt_specific", "Cic_specific",
           title_str = "TEST: Runt-specific vs Cic-specific DARs",
           fname = "TEST_venn_runt_vs_cic.pdf",
           out_dir = output_dir)
cat("[PASS] Venn diagram written to ", file.path(output_dir, "TEST_venn_runt_vs_cic.pdf"), "\n")

cat("\n================ ALL TESTS PASSED ================\n")
cat("Check the following files in", output_dir, "to visually confirm output:\n")
cat("  - Explanatory_Runt_specific_panel.pdf\n")
cat("  - Explanatory_Cic_specific_panel.pdf\n")
cat("  - TEST_venn_runt_vs_cic.pdf  (should show two circles overlapping by 10)\n")
