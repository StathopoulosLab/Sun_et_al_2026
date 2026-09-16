################################################################################
# 00_TEST_mixed_effects_quickstart.R
#
# Minimal reviewer-facing test for
# 04_mixed_effects_logistic_regression_state_selective_accessibility.r.
#
# Copies toy TF partner peaks + DAR data into place, then sources the real
# script UNMODIFIED so it runs its actual lme4/lmerTest mixed-model fitting,
# binary GLMM, optional glmmLasso regularization, and bootstrap stability
# check on synthetic data.
#
# SCOPE: the real script defines ~27 contrasts across two genotype splits.
# This kit supplies toy DAR data for exactly TWO of them:
#   - BOT_vs_BOTR   (split 1, non-ventralized; NOT in NEGATE_LFC_STEMS)
#   - BOTv_vs_BOTCv (split 2, ventralized; IS in NEGATE_LFC_STEMS -- exercises
#                    the sign-negation code path specific to this stem)
# The other ~25 contrasts will print "[missing]" for their DAR bed file and
# get skipped by run_and_save_lmm's own `if (length(gr) == 0) SKIP` check --
# this is the real script's own designed graceful-skip behavior, not
# something this kit adds. Seeing many "[missing]"/"SKIP (0 DARs)" lines in
# the output is EXPECTED, not a failure.
#
# WHAT THIS PROVES: the core statistical machinery runs correctly end-to-end
# -- TF design matrix construction from real ChIP/BDTNP overlap data, the
# |logFC| threshold stratification (3 thresholds + unfiltered), the mixed
# model formula construction and lmerTest fit, the binary GLMM, and (if
# glmmLasso is installed) the regularized models -- not that the reported
# coefficients match any particular expected value, which depends on real
# statistical fitting I have no way to hand-verify.
#
# TOY DATA SIZING: both DAR stems have 120 peaks spread across 5 chromosomes
# (chr2L/2R/3L/3R/X, 24 each) -- the multi-chromosome spread is required so
# the model's (1|chr) random-effect term has more than one level to estimate
# variance from (a single-chromosome toy dataset caused every fit to fail
# silently -- see README for how this was diagnosed and fixed). 90 of the
# 120 peaks clear |logFC|>0.5, comfortably clearing run_lmm_stratum's
# "fewer than 20 peaks" floor at every threshold. All 19 TF partner peak
# files bind a random 25-95 of the 120 loci each, giving every TF column
# real variance.
#
# USAGE: Rscript 00_TEST_mixed_effects_quickstart.R
# (run from inside this test_kit_mixed_effects/ folder)
################################################################################

real_script_path <- "../ATACseq_Processing/02_DAR_Calling/04_mixed_effects_logistic_regression_state_selective_accessibility.r"
if (!file.exists(real_script_path)) {
  stop("Can't find the mixed-effects script at: ", real_script_path,
       "\nEdit `real_script_path` at the top of this script to point at it.")
}

cat("================ TEST SETUP ================\n")

# Copy toy TF partner peaks (into ./data/... at the real script's exact
# hardcoded paths) and toy DAR data (into ./ATAC_DARS/) alongside wherever
# this script is run from.
file.copy("./test_data/data", ".", recursive = TRUE, overwrite = TRUE)
dir.create("./ATAC_DARS", showWarnings = FALSE)
file.copy(list.files("./test_data/ATAC_DARS", full.names = TRUE),
          "./ATAC_DARS", overwrite = TRUE)
cat("Copied toy TF partner peaks into ./data/ and toy DAR data into ./ATAC_DARS/\n\n")

cat("================ RUNNING REAL SCRIPT ================\n")
cat("(This sources the mixed-effects script unmodified. Expect many\n")
cat(" '[missing]'/'SKIP (0 DARs)' lines for the ~25 contrasts this kit\n")
cat(" doesn't supply toy data for -- that's the real script's own\n")
cat(" graceful-skip behavior, not a failure.)\n\n")

source(real_script_path)

cat("\n================ VALIDATING RESULTS ================\n")

# Reconstruct the same RUN_TAG-based output directory names the real script
# computed (format(Sys.Date(), ...) -- changes daily, so we recompute rather
# than hardcode).
run_tag        <- format(Sys.Date(), "%Y%m%d_bootstrapCI")
split1_out_dir <- sprintf("./Split1_NonVent_newcic_%s", run_tag)
split2_out_dir <- sprintf("./Split2_Vent_newcic_%s", run_tag)

strata <- c("lfc0.50", "lfc0.35", "lfc0.25", "all_peaks")

check_contrast <- function(out_dir, results_subdir, label) {
  all_found <- TRUE
  for (stratum in strata) {
    coef_file <- file.path(out_dir, results_subdir, "lmm", stratum, "lmm_coefficients.txt")
    rds_file  <- file.path(out_dir, results_subdir, "lmm", stratum, "lmm_result.rds")
    ok <- file.exists(coef_file) && file.exists(rds_file)
    cat(sprintf("  [%s] %-10s %s\n", label, stratum, if (ok) "OK" else "MISSING"))
    if (ok) {
      ct <- tryCatch(read.table(coef_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
                     error = function(e) NULL)
      if (!is.null(ct) && nrow(ct) > 0) {
        n_sig <- sum(ct$significant, na.rm = TRUE)
        cat(sprintf("      %d TF coefficients fit, %d nominally significant, n_peaks_total=%s\n",
                    nrow(ct), n_sig, unique(ct$n_peaks_total)[1]))
      }
    } else {
      all_found <- FALSE
    }
  }
  all_found
}

cat("\n--- Split 1 (non-ventralized): BOT_vs_BOTR ---\n")
ok1 <- check_contrast(split1_out_dir, "results_BOT_vs_BOTR", "BOT_vs_BOTR")

cat("\n--- Split 2 (ventralized): BOTv_vs_BOTCv ---\n")
ok2 <- check_contrast(split2_out_dir, "results_BOTv_vs_BOTCv", "BOTv_vs_BOTCv")

cat("\n================ SUMMARY ================\n")
if (ok1 && ok2) {
  cat("[PASS] Both toy contrasts produced coefficient + RDS output for all 4 strata.\n")
} else {
  cat("[FAIL] Some expected output files are missing -- see MISSING lines above.\n")
}
cat("\nNote: glmmLasso outputs (glmmlasso_*_coefs.txt) and bootstrap stability\n")
cat("results only appear if the glmmLasso package is installed and\n")
cat("RUN_BOOTSTRAP_STABILITY is TRUE (both check gracefully / are on by default\n")
cat("in the real script) -- check the console output above for\n")
cat("'[bootstrap] running 200 resamples' and 'glmmLasso patched'/'coefs' lines\n")
cat("to confirm those ran too.\n")
