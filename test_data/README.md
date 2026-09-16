# Test data

Small synthetic datasets that let a reviewer run key parts of the analysis
pipeline without needing the real (large, restricted) sequencing data. Each
subfolder is self-contained: a README explaining what's simulated and what
output to expect, a script to generate/copy the test data into place, and a
script to run the real pipeline code against it.

| Folder | Tests | Run with |
|---|---|---|
| `test_kit_tf_enrichment/` | `03_TF_Enrichment/00_tf_venn_shared_helpers.R` | `Rscript 00_TEST_tf_enrichment_quickstart.R` |
| `test_kit_dar_calling/` | `02_DAR_Calling/01_dar_calling_seven_genotype_nc14b.r` | `Rscript generate_toy_dar_cache.R` then the real script |
| `test_kit_peak_fate/` | `05_Peak_Fate/01_peak_fate_classification.r` | `Rscript 00_TEST_peak_fate_quickstart.R` |
| `test_kit_rnaseq_deg/` | `RNAseq_Processing/Analysis_DEG_LimmaBOTvBOTCv.R` | `Rscript 00_TEST_rnaseq_deg_quickstart.R` |

None of these modify the actual analysis scripts — each kit only supplies
synthetic input data in the same format the real pipeline expects. See the
README inside each subfolder for what's simulated, expected output, and any
known limitations of that specific test.
