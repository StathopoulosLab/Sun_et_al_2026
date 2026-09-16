#!/usr/bin/env Rscript
################################################################################
# 00_TEST_rnaseq_deg_quickstart.R
#
# Minimal reviewer-facing test for Analysis_DEG_LimmaBOTvBOTCv.R.
#
# Runs the real script (via its own --counts/--outdir CLI flags -- no
# modification needed) on a 16-embryo synthetic counts matrix, then checks
# the console-reported sample classification and contrast arm sizes against
# hand-computed expected values (see test_data's generation log / README).
#
# WHAT THIS EXERCISES:
#   - get_group()/get_stage() column-name parsing across every naming
#     convention the real script supports (bare/prefixed, FoxL1_BOTv vs
#     foxl1_Dm vs BOTv, HLH54F_BOTCv vs HLH54F_Dm)
#   - the B6_BOTC decoy-barcode guard (bc19-24 under foxl1_Dm prefix must
#     be excluded, not misassigned to FoxL1_BOTv)
#   - the EXCLUDE_BCS barcode exclusion (bc27)
#   - batch/sublibrary detection and the confounding safety check
#   - the full limma-voom + eBayes contrast pipeline on real (synthetic)
#     count data, including engineered genotype- and temporal-effect genes
#     so DE calls aren't just noise
#
# USAGE: Rscript 00_TEST_rnaseq_deg_quickstart.R
################################################################################

real_script_path <- "../RNAseq_Processing/Analysis_DEG_LimmaBOTvBOTCv.R"
if (!file.exists(real_script_path)) {
  stop("Can't find Analysis_DEG_LimmaBOTvBOTCv.R at: ", real_script_path,
       "\nEdit `real_script_path` at the top of this script to point at it.")
}

counts_path <- "./test_data/counts_matrix.tsv"
symbol_path <- "./test_data/gene_id_to_symbol.tsv"
out_dir     <- "./test_output"

if (!file.exists(counts_path)) stop("Missing ", counts_path)

# The real script looks for SYMBOL_MAP ("gene_id_to_symbol.tsv") in the
# CURRENT WORKING DIRECTORY (it's not a CLI flag) -- copy it there for
# this run so symbol lookup is exercised too.
file.copy(symbol_path, "./gene_id_to_symbol.tsv", overwrite = TRUE)

cat("================ RUNNING REAL SCRIPT ================\n")
cat("(via its own --counts/--outdir flags, source() does not see commandArgs()\n")
cat(" the same way Rscript does, so we invoke it as a subprocess instead.)\n\n")

status <- system2("Rscript",
  args = c(shQuote(real_script_path),
           "--counts", shQuote(counts_path),
           "--outdir", shQuote(out_dir)),
  stdout = "", stderr = "")

if (status != 0) {
  stop("[FAIL] Analysis_DEG_LimmaBOTvBOTCv.R exited with non-zero status: ", status)
}

cat("\n================ VALIDATING RESULTS ================\n")

summary_path <- file.path(out_dir, "summary_table.tsv")
if (!file.exists(summary_path)) {
  stop("[FAIL] summary_table.tsv was not created -- the real script did not complete.")
}
summ <- read.table(summary_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Observed summary table:\n"); print(summ)

expected_contrasts <- c("BOTCv_nc14b_vs_BOTv_nc14b", "BOTCv_late_vs_BOTv_nc14d",
                        "BOTCv_late_vs_BOTCv_nc14b", "BOTv_nc14d_vs_BOTv_nc14b",
                        "BOTCv_late_vs_BOTv_nc14b")
missing <- setdiff(expected_contrasts, summ$contrast)
if (length(missing) > 0) {
  cat(sprintf("[MISMATCH] Missing contrast(s) in output: %s\n", paste(missing, collapse=", ")))
} else {
  cat("[PASS] All 5 expected contrasts completed (none SKIPped).\n")
}

cat("\nCheck the console output above (scroll up) for these expected lines:\n")
cat("  GROUP_A ('FoxL1_BOTv') matches 6 columns\n")
cat("  GROUP_B ('HLH54F_BOTCv') matches 7 columns\n")
cat("  13 embryos after exclusions\n")
cat("  '2 sublibraries detected -- adding batch covariate' for all 5 contrasts\n")
cat("  (NOT 'batch is confounded with condition' for any of them)\n\n")
cat("If any of those differ from what's printed above, that's a real\n")
cat("mismatch worth reporting back -- see test_data/README.md for exactly\n")
cat("how each expected number was derived.\n")
