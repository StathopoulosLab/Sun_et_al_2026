# Toy test data for TF-enrichment helper functions

This is **synthetic data**, not real ChIP-seq/ATAC-seq peaks. It exists only
to let a reviewer confirm that `00_tf_venn_shared_helpers.R`'s functions run
correctly and produce sensible output, without needing the real (large,
restricted) sequencing data.

## How it was built

All peaks are 200bp, on `chr2L`, spaced 3000bp apart starting at position
1,000,000, indexed 0, 1, 2, ... Two peaks with the same index have **identical
coordinates** (guaranteed overlap); different indices never overlap.

| Set | Peak indices used | Count |
|---|---|---|
| `DAR_Runt_specific.bed` | 0–9 (shared) + 10–29 (unique) | 30 |
| `DAR_Cic_specific.bed`  | 0–9 (shared) + 30–44 (unique) | 25 |
| `TF_Runt.bed`  | 0–14 | 15 |
| `TF_Cic.bed`   | 30–41 | 12 |
| `TF_Zld.bed`   | 5–9, 20–24, 35–39 | 15 |
| `TF_Bcd.bed`   | 50–54 (negative control — doesn't overlap either DAR set) | 5 |

## Expected results when run through the real pipeline functions

**`tf_explanatory_summary()` on Runt_specific (n=30):**
- Explained by ≥1 TF: 20 (66.7%)
- Single-TF explained: 15
- Multi-TF explained: 5
- Unexplained: 10

**`tf_explanatory_summary()` on Cic_specific (n=25):**
- Explained by ≥1 TF: 22 (88.0%)
- Single-TF explained: 12
- Multi-TF explained: 10
- Unexplained: 3

**`make_venn2()` on Runt_specific vs Cic_specific:**
- Shared peaks (genomic overlap): 10 — matches the 10 peaks both sets were
  built to share

These numbers are asserted automatically by `stopifnot()` calls in
`00_TEST_tf_enrichment_quickstart.R` — if the script prints `ALL TESTS PASSED`,
the pipeline logic is verified end-to-end. `TF_Bcd` is included specifically
as a negative control: it should contribute 0 to any explained-DAR count for
either set, since none of its 5 peaks fall inside the DAR test regions.
