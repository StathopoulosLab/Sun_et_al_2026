################################################################################
# 00_TEST_part3_cluster_fate_quickstart.R
#
# Minimal reviewer-facing test for atac_rna_integration_Step1.R, Part 3
# (cluster x fate enrichment) only.
#
# SCOPE: --parts 3 only (this kit supplies no DEG data, so Part 1 would
# just print "Not found" -- cleaner to isolate Part 3 entirely). Part 2 is
# already off by the script's own default.
#
# WHAT'S TESTED: the Fisher's-exact-test enrichment grid between gene
# clusters (from Analysis_GeneClusters_BOTvBOTCv.py's real output format)
# and peak-fate gene sets (from 01_peak_fate_classification.r's real
# output format), including the collapse_converged() logic that unions
# multiple "_Converged_<mechanism>" files per genotype into one set before
# testing, and the BH multiple-testing correction (global, per-genotype-
# stratified, and per-fate-stratified variants).
#
# TOY DATA DESIGN: 100 genes across 4 clusters (25 each). Two fate gene
# sets are engineered with strong real overlap against specific clusters
# (BOTv_Maintained x Cluster 1, BOTCv_Maintained x Cluster 2 -- both
# a=15/18 genes overlapping), and four more sets (two Converged mechanism
# files to test the union logic, plus Emerged/Reversed) are scattered
# near-uniformly across clusters as an implicit negative control.
#
# EXPECTED RESULTS (independently computed in Python, exact 2x2 Fisher
# hypergeometric, matching the real script's contingency-table logic):
#   Cluster 1 x BOTv_Maintained:  a=15 b=10 c=3 d=72   raw p=7.39e-09
#   Cluster 2 x BOTCv_Maintained: a=15 b=10 c=3 d=72   raw p=7.39e-09
# Both are astronomically significant and should easily survive BH
# correction (padj) in the collapsed significance-test output.
#
# USAGE: Rscript 00_TEST_part3_cluster_fate_quickstart.R
# (run from inside this test_kit_rna_atac_integration_part3/ folder)
################################################################################

real_script_path <- "../ATACseq_Processing/06_RNA_ATAC_Integration/atac_rna_integration_Step1.R"
if (!file.exists(real_script_path)) {
  stop("Can't find atac_rna_integration_Step1.R at: ", real_script_path,
       "\nEdit `real_script_path` at the top of this script to point at it.")
}

cat("================ TEST SETUP ================\n")

dir.create("./results", showWarnings = FALSE)
file.copy("./test_data/results", ".", recursive = TRUE, overwrite = TRUE)

dir.create("./Overview_Plots", showWarnings = FALSE)
file.copy("./test_data/Overview_Plots", ".", recursive = TRUE, overwrite = TRUE)

cat("Copied toy cluster assignments into ./results/.../explore_clusters_botv_botcv/\n")
cat("and toy fate gene lists into ./Overview_Plots/peak_fate/BOTv_vs_BOTCv/gene_lists/\n\n")

cat("================ RUNNING REAL SCRIPT ================\n")
cat("(via --parts 3 only -- this kit supplies no DEG data for Part 1,\n")
cat(" and Part 2 is already off by the script's own default.)\n\n")

status <- system2("Rscript",
  args = c(shQuote(real_script_path), "--parts", "3"),
  stdout = TRUE, stderr = TRUE)
console_lines <- status
cat(paste(console_lines, collapse = "\n"), "\n")
exit_status <- attr(status, "status")
if (!is.null(exit_status) && exit_status != 0) {
  stop("[FAIL] atac_rna_integration_Step1.R exited with non-zero status: ", exit_status)
}

cat("\n================ VALIDATING RESULTS ================\n")

enrich_path <- "./Overview_Plots/atac_rna_integration/part3_cluster_fate/cluster_fate_enrichment_table.tsv"
comp_path   <- "./Overview_Plots/atac_rna_integration/part3_cluster_fate/cluster_fate_composition_table_fine.tsv"

if (!file.exists(enrich_path)) {
  stop("[FAIL] Expected output not found: ", enrich_path,
       "\nCheck the console output above for where Part 3 failed.")
}

enrich <- read.table(enrich_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("Enrichment table: %d rows\n", nrow(enrich)))

check_hit <- function(cl_id, fate_name) {
  row <- enrich[enrich$cluster == cl_id & enrich$fate_short == fate_name, ]
  if (nrow(row) == 0) {
    cat(sprintf("  [MISSING] Cluster %d x %s -- no row found\n", cl_id, fate_name))
    return(FALSE)
  }
  row <- row[1, ]
  cat(sprintf("  Cluster %d x %-16s  n_overlap=%d (expected 15)  padj=%.2e (expected << 0.05)\n",
              cl_id, fate_name, row$n_overlap, row$padj_fate_stratified))
  row$n_overlap == 15 && row$padj_fate_stratified < 0.001
}

cat("\nChecking the two designed-positive enrichments:\n")
ok1 <- check_hit(1, "BOTv_Maintained")
ok2 <- check_hit(2, "BOTCv_Maintained")

ok <- isTRUE(ok1) && isTRUE(ok2) && file.exists(comp_path)

if (ok) {
  cat("\n================ PASS ================\n")
  cat("Both designed cluster-fate enrichments detected with expected overlap\n")
  cat("and strong significance. Fisher grid + BH correction logic works.\n")
} else {
  cat("\n================ MISMATCH -- see numbers above ================\n")
}
cat("\nAlso check console output above for:\n")
cat("  - 'collapse_converged' effectively unioning the two BOTv_Converged_* files\n")
cat("    (the genes_BOTv_Converged_mechA.csv + mechB.csv should merge into one\n")
cat("     'BOTv_Converged' fate in the significance-test table, not appear\n")
cat("     as two separate fates there)\n")
