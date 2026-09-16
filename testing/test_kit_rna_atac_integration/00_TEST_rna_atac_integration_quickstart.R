################################################################################
# 00_TEST_rna_atac_integration_quickstart.R
#
# Minimal reviewer-facing test for atac_rna_integration_Step1.R, Part 1
# (concordance bins) only.
#
# SCOPE: --parts 1 only. Part 2 (GSEA-style enrichment) is already excluded
# from the real script's own default (--parts defaults to "1,3", and this
# kit narrows that further to just "1"). Part 3 (cluster x fate enrichment)
# needs CLUSTER_DIR, which still points at an external path outside this
# repo (see the script's own TODO comment in CONFIG) -- not covered by this
# kit. This kit exercises the core ATAC x RNA merge/concordance logic that
# both Part 1 and Part 3 share, without needing cluster-assignment data.
#
# WHAT'S TESTED: exactly 1 of the 5 RNA_CONTRASTS x their ATAC pairs
# (BOTCv_nc14b_vs_BOTv_nc14b x BOTCv_vs_BOTv_nc14b, i.e. stem BOTv_vs_BOTCv).
# The script's own driver loop gracefully `next`s past RNA contrasts or ATAC
# stems with no data file (read_rna()/read_atac_annotated() return NULL) --
# so running the real script as-is against this kit's toy data will process
# the 1 contrast we supply and silently skip the other 4. That's the real
# script's own designed behavior, not something this kit adds.
#
# TOY DATA DESIGN: 60 genes with IDENTICAL symbols in both the toy RNA
# results_full.tsv and the toy ATAC DAR annotated file, guaranteeing
# merge_atac_rna()'s symbol-based merge succeeds (real biological data
# would need ChIPseeker/gene-annotation packages for this; toy data sidesteps
# that entirely by construction). Engineered as:
#   - 25 genes: concordant (ATAC and RNA same direction, both clear their
#     respective significance/effect-size thresholds)
#   - 15 genes: discordant (opposite directions, both clear thresholds)
#   - 20 genes: weak ATAC signal (|log2FC| < 0.50, won't clear LFC_ATAC
#     threshold, exercises the "faded" non-significant point rendering)
#
# EXPECTED RESULTS:
#   40/60 RNA genes have padj<0.10 (sig_rna)
#   40/60 ATAC peaks have |log2FC|>0.50 (sig_atac)
#   All 60 genes merge successfully (n=60 in the merged concordance table)
#   Roughly 25 Concordant, ~15+20 split across Discordant/other quadrants
#     depending on exact random sign draws -- see test_data's generation
#     log for the precise designed split
#
# USAGE: Rscript 00_TEST_rna_atac_integration_quickstart.R
# (run from inside this test_kit_rna_atac_integration/ folder)
################################################################################

real_script_path <- "../ATACseq_Processing/06_RNA_ATAC_Integration/atac_rna_integration_Step1.R"
if (!file.exists(real_script_path)) {
  stop("Can't find atac_rna_integration_Step1.R at: ", real_script_path,
       "\nEdit `real_script_path` at the top of this script to point at it.")
}

cat("================ TEST SETUP ================\n")

dir.create("./ATAC_DARS", showWarnings = FALSE)
file.copy(list.files("./test_data/ATAC_DARS", full.names = TRUE),
          "./ATAC_DARS", overwrite = TRUE)

dir.create("./results", showWarnings = FALSE)
file.copy("./test_data/results", ".", recursive = TRUE, overwrite = TRUE)

cat("Copied toy DAR data into ./ATAC_DARS/ and toy DEG output into ./results/\n")
cat("(matches DEG_DIR's fixed default: results/combined_figures/qc_and_deg_botv_botcv/deg_limma)\n\n")

cat("================ RUNNING REAL SCRIPT ================\n")
cat("(via --parts 1, skipping Part 2 [already off by default] and Part 3\n")
cat(" [needs CLUSTER_DIR, not covered by this kit]. Expect 'Not found' /\n")
cat(" skip messages for the other 4 RNA_CONTRASTS this kit doesn't supply\n")
cat(" data for -- that's the real script's own graceful-skip behavior.)\n\n")

status <- system2("Rscript",
  args = c(shQuote(real_script_path), "--parts", "1"),
  stdout = TRUE, stderr = TRUE)
console_lines <- status
cat(paste(console_lines, collapse = "\n"), "\n")
exit_status <- attr(status, "status")
if (!is.null(exit_status) && exit_status != 0) {
  stop("[FAIL] atac_rna_integration_Step1.R exited with non-zero status: ", exit_status)
}

cat("\n================ VALIDATING RESULTS ================\n")

out_stem <- "./Overview_Plots/atac_rna_integration/part1_concordance/concordance_BOTCv_nc14b_vs_BOTv_nc14b_BOTCv_vs_BOTv_nc14b"
tsv_path <- paste0(out_stem, ".tsv")
pdf_path <- paste0(out_stem, ".pdf")

# NOTE: we deliberately do NOT read the .tsv back with read.table() here.
# The real script writes its quadrant labels with a literal embedded "\n"
# (e.g. "Concordant\nOpen+Up", used for legend line-breaking in the plot),
# via write.table(..., quote=FALSE) -- so that newline lands as a raw line
# break in the file, splitting one data row across two physical lines.
# Any standard TSV parser chokes on this (confirmed in an actual test run).
# Instead we parse the real script's own console output, which already
# reports the exact same numbers in an unambiguous form.
grep_num <- function(pattern) {
  hit <- grep(pattern, console_lines, value = TRUE)
  if (length(hit) == 0) return(NA_integer_)
  m <- regmatches(hit[1], regexec(pattern, hit[1]))[[1]]
  if (length(m) < 2) return(NA_integer_)
  as.integer(m[2])
}

n_merged      <- grep_num("Merged: ([0-9]+) genes")
n_conc_close  <- grep_num("Concordant Close\\+Down: ([0-9]+)")
n_conc_open   <- grep_num("Concordant Open\\+Up: ([0-9]+)")
n_disc_close  <- grep_num("Discordant Close\\+Up: ([0-9]+)")
n_disc_open   <- grep_num("Discordant Open\\+Down: ([0-9]+)")

cat(sprintf("Merged table: %s genes (expected: 60)\n", n_merged))
cat("\nQuadrant breakdown (expected: 19/17/12/12):\n")
cat(sprintf("  Concordant Open+Up:    %s\n", n_conc_open))
cat(sprintf("  Concordant Close+Down: %s\n", n_conc_close))
cat(sprintf("  Discordant Open+Down:  %s\n", n_disc_open))
cat(sprintf("  Discordant Close+Up:   %s\n", n_disc_close))

ok <- isTRUE(!is.na(n_merged) && n_merged == 60 &&
  identical(c(n_conc_open, n_conc_close, n_disc_open, n_disc_close), c(19L, 17L, 12L, 12L)) &&
  file.exists(pdf_path) && file.exists(tsv_path))

if (ok) {
  cat("\n================ PASS ================\n")
  cat("60 genes merged, exact expected quadrant split, plot + table saved.\n")
  cat("Core ATAC x RNA concordance logic works.\n")
} else {
  cat("\n================ MISMATCH -- see numbers above ================\n")
}
cat("\nVisual check: ", pdf_path, "\n")
cat("(should show a scatter of ATAC log2FC vs RNA log2FC, colored by quadrant,\n")
cat(" plus a quadrant-count bar and an annotation-feature breakdown panel)\n")
