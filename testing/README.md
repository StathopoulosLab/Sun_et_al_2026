# Test data

Small synthetic datasets that let a reviewer run key parts of the analysis
pipeline without needing the real (large, restricted) sequencing data. Each
subfolder is self-contained: a README explaining what's simulated and what
output to expect, a script to generate/copy the test data into place, and a
script to run the real pipeline code against it.

All seven kits below have been run against the real scripts and confirmed
passing. Three real bugs were found and fixed along the way (see each
kit's own README for details): a hyphen/dot column-name mismatch in the
DAR-calling cache, a degenerate-random-effect issue in the mixed-effects
toy data, and a ggplot2 size-legend crash in the RNA/ATAC integration
script's cluster-fate dotplot.

| Folder | Tests | Run with |
|---|---|---|
| `test_kit_tf_enrichment/` | `03_TF_Enrichment/00_tf_venn_shared_helpers.R` | `Rscript 00_TEST_tf_enrichment_quickstart.R` |
| `test_kit_dar_calling/` | `02_DAR_Calling/01_dar_calling_seven_genotype_nc14b.r` | `Rscript generate_toy_dar_cache.R` then the real script |
| `test_kit_peak_fate/` | `05_Peak_Fate/01_peak_fate_classification.r` | `Rscript 00_TEST_peak_fate_quickstart.R` |
| `test_kit_rnaseq_deg/` | `RNAseq_Processing/Analysis_DEG_LimmaBOTvBOTCv.R` | `Rscript 00_TEST_rnaseq_deg_quickstart.R` |
| `test_kit_mixed_effects/` | `02_DAR_Calling/04_mixed_effects_logistic_regression_state_selective_accessibility.r` | `Rscript 00_TEST_mixed_effects_quickstart.R` |
| `test_kit_rna_atac_integration/` | `06_RNA_ATAC_Integration/atac_rna_integration_Step1.R` (Part 1: concordance) | `Rscript 00_TEST_rna_atac_integration_quickstart.R` |
| `test_kit_rna_atac_integration_part3/` | `06_RNA_ATAC_Integration/atac_rna_integration_Step1.R` (Part 3: cluster × fate enrichment) | `Rscript 00_TEST_part3_cluster_fate_quickstart.R` |

None of these modify the actual analysis scripts — each kit only supplies
synthetic input data in the same format the real pipeline expects. See the
README inside each subfolder for what's simulated, expected output, and any
known limitations of that specific test.

Not covered (lower priority — see individual kit READMEs and prior review
notes for rationale): the other two DAR-calling script variants, the
strict-mode peak-fate script, Part 2 of the integration script (GSEA-style
rank enrichment, off by the real script's own default), and plotting-only
or HPC-preprocessing scripts generally.
