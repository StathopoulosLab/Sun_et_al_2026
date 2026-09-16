################################################################################
# 00_TEST_peak_fate_quickstart.R
#
# Minimal reviewer-facing test for 01_peak_fate_classification.r.
#
# Copies the toy DAR files into place, then sources the real script so it
# runs its ACTUAL fate-classification logic (Part A: Maintained/Reversed/
# Converged/Emerged; Part B: non-vent Gain/Lose) on synthetic data with
# known correct answers.
#
# WHAT'S DELIBERATELY OMITTED (both are fully optional -- confirmed by
# tracing the script's guard conditions):
#   - called-peak support files (Merged_Peak_BEDs/) -- Emerged-locus
#     rescue/exclusion tiers are skipped; script prints [NOTE] and continues
#   - *_FULL_results.txt unfiltered result files -- trend-rescue fallback
#     paths for Converged/Emerged are skipped; script prints [NOTE] and continues
#   - ChIPseeker / TxDb.Dmelanogaster / org.Dm.eg.db -- if not installed,
#     annotation plots/gene-window lookups are skipped entirely (has_chipseeker
#     check at the top of the real script)
# None of these omissions should cause an error -- if one does, that's a bug
# worth reporting, not expected behavior.
#
# EXPECTED RESULTS (Part A, BOTv vs BOTCv):
#   Maintained = 8, Reversed = 4, Converged = 8, Emerged = 5
#   fate_df should have 25 rows total (20 nc14b-based + 5 Emerged)
#
# EXPECTED RESULTS (Part B, non-vent, each of 3 genotype pairs):
#   10 temporal DARs: 5 gain / 5 lose accessibility
#   5 of the 10 also overlap that genotype's nc14b anchor DAR set (also_nc14b=TRUE)
#
# USAGE: Rscript 00_TEST_peak_fate_quickstart.R
# (run from inside this test_kit/ folder)
################################################################################

real_script_path <- "../ATACseq_Processing/05_Peak_Fate/01_peak_fate_classification.r"
if (!file.exists(real_script_path)) {
  stop("Can't find 01_peak_fate_classification.r at: ", real_script_path,
       "\nEdit `real_script_path` at the top of this script to point at it.")
}

cat("================ TEST SETUP ================\n")
dir.create("./ATAC_DARS", showWarnings = FALSE)
file.copy(list.files("./test_data/ATAC_DARS", full.names = TRUE),
          "./ATAC_DARS", overwrite = TRUE)
cat("Copied toy DAR files into ./ATAC_DARS/\n\n")

cat("================ RUNNING REAL SCRIPT ================\n")
cat("(This sources 01_peak_fate_classification.r unmodified -- its own\n")
cat(" config points at ./ATAC_DARS, ./Overview_Plots relative to THIS directory.)\n\n")

# The real script's CONFIG section hardcodes dar_dir <- "./ATAC_DARS" and
# out_dir <- "./Overview_Plots" -- both relative to the working directory,
# which is this test_kit folder when run this way. No modification needed.
source(real_script_path)

cat("\n================ VALIDATING RESULTS ================\n")

fate_csv <- file.path("./Overview_Plots/peak_fate/BOTv_vs_BOTCv", "peak_fate_data.csv")
if (!file.exists(fate_csv)) {
  stop("[FAIL] peak_fate_data.csv was not created -- the real script did not complete Part A.")
}
fate_out <- read.csv(fate_csv, stringsAsFactors = FALSE)
observed <- table(fate_out$fate)
cat("Observed fate counts:\n"); print(observed)

expected <- c(Maintained = 8, Reversed = 4, Converged = 8, Emerged = 5)
cat("\nExpected fate counts:\n"); print(expected)

mismatch <- FALSE
for (f in names(expected)) {
  obs_n <- if (f %in% names(observed)) observed[[f]] else 0
  if (obs_n != expected[[f]]) {
    cat(sprintf("[MISMATCH] %s: expected %d, got %d\n", f, expected[[f]], obs_n))
    mismatch <- TRUE
  }
}

nv_csv <- file.path("./Overview_Plots/peak_fate/NonVent", "nonvent_temporal_data.csv")
if (file.exists(nv_csv)) {
  nv_out <- read.csv(nv_csv, stringsAsFactors = FALSE)
  cat("\nNon-vent rows written:", nrow(nv_out), " (expected: 30 = 10 x 3 genotype pairs)\n")
  if (nrow(nv_out) != 30) mismatch <- TRUE
} else {
  cat("\n[FAIL] nonvent_temporal_data.csv was not created -- Part B did not complete.\n")
  mismatch <- TRUE
}

if (mismatch) {
  cat("\n================ SOME CHECKS FAILED — see [MISMATCH]/[FAIL] above ================\n")
} else {
  cat("\n================ ALL CHECKS PASSED ================\n")
}
cat("\nInspect these for a visual check:\n")
cat("  ./Overview_Plots/peak_fate/BOTv_vs_BOTCv/peak_fate_alluvial.pdf (if ggalluvial installed)\n")
cat("  ./Overview_Plots/peak_fate/BOTv_vs_BOTCv/peak_fate_stackedbar.pdf\n")
cat("  ./Overview_Plots/peak_fate/NonVent/nonvent_temporal_fate.pdf\n")
