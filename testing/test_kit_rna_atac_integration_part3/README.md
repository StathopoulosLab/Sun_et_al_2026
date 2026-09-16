# Test kit: RNA/ATAC integration Part 3 (cluster × fate enrichment)

Runs your **actual, unmodified** Part 3 logic — Fisher's-exact-test
enrichment between gene expression clusters and ATAC peak-fate categories,
including the `collapse_converged()` mechanism-union step and three
flavors of BH multiple-testing correction — on synthetic data with a
precisely known correct answer.

## Prerequisite fix (already applied, separately)

`CLUSTER_DIR` previously pointed at an external path
(`../../seRNA-seq_03_19_26/...`) outside this repo, the same problem
`DEG_DIR` had. Fixed to match `Analysis_GeneClusters_BOTvBOTCv.py`'s own
default `--outdir`. This kit assumes that fix is already in place.

## Scope

This kit runs **only Part 3** (`--parts 3`) — it supplies no DEG data, so
Part 1 would just print `"Not found"` noise; cleaner to isolate Part 3
entirely. Part 2 is already off by the script's own default.

## What's simulated

**`cluster_assignments.tsv`** (matching `Analysis_GeneClusters_BOTvBOTCv.py`'s
real 5-column output): 100 genes evenly split across 4 clusters (25 each) —
`Early-responsive`, `Late-responsive`, `Stable-high`, `Stable-low`.

**Fate gene lists** (matching `01_peak_fate_classification.r`'s real
`genes_<fate>.csv` naming and `SYMBOL` column format), 6 files:

| File | Genes | Design |
|---|---|---|
| `genes_BOTv_Maintained.csv` | 18 | **15 from Cluster 1** — real designed enrichment |
| `genes_BOTCv_Maintained.csv` | 18 | **15 from Cluster 2** — real designed enrichment |
| `genes_BOTv_Converged_mechA.csv` | 10 | Scattered — tests `collapse_converged()` union |
| `genes_BOTv_Converged_mechB.csv` | 8 | Scattered — unions with mechA into one `BOTv_Converged` set |
| `genes_BOTCv_Emerged.csv` | 12 | Scattered — implicit negative control |
| `genes_BOTv_Reversed.csv` | 10 | Scattered — implicit negative control |

## Expected results

Computed independently in Python using the exact same 2×2 contingency
table logic the real script's `run_fisher_grid()` uses (hypergeometric
right-tail test):

```
Cluster 1 x BOTv_Maintained:  a=15 b=10 c=3 d=72   raw p=7.39e-09
Cluster 2 x BOTCv_Maintained: a=15 b=10 c=3 d=72   raw p=7.39e-09
```

Both are astronomically significant and should easily survive BH
correction — the test script checks for `n_overlap=15` and
`padj_fate_stratified < 0.001` on both.

Also worth checking in the console output: the two `BOTv_Converged_mechA`/
`mechB` files should merge into a single `BOTv_Converged` entry in the
significance-test table (not appear as two separate fates there) —
confirms `collapse_converged()`'s union logic is working, including the
`"genes_"`-prefix-stripping fix documented in the real script's own
comments (a past bug where that prefix wasn't stripped before matching,
silently disabling the collapse entirely).

## How to install and run

1. Copy this whole `test_kit_rna_atac_integration_part3/` folder to the
   root of `Sun_et_al_2026/` (sibling of `ATACseq_Processing/`).
2. From inside `test_kit_rna_atac_integration_part3/`, run:
   ```
   Rscript 00_TEST_part3_cluster_fate_quickstart.R
   ```

Requires: same packages as the Part 1 kit (`GenomicRanges`, `ggplot2`,
`dplyr`, `patchwork`, `scales`); `tidyr` (optional, used by some Part 3
sub-plots — degrades gracefully if absent, per the script's own
`has_tidyr` check).

## What success looks like

```
================ PASS ================
Both designed cluster-fate enrichments detected with expected overlap
and strong significance. Fisher grid + BH correction logic works.
```

Output lands in `./Overview_Plots/atac_rna_integration/part3_cluster_fate/`
— the enrichment and composition TSVs plus whatever plots Part 3
generates downstream of the Fisher grid.

## Important: I have not executed this script

Same caveat as every other kit — no R available where I built this. I
traced `run_fisher_grid()`'s exact contingency-table construction and
independently verified the p-values in Python, and traced
`collapse_converged()`'s union logic including the prefix-stripping fix
noted in the real script's own comments — but Part 3 is long (the
Fisher-grid section alone is a few hundred lines, with plotting code
beyond what's shown here that I haven't fully traced). **Please run this
and paste the full console output.**

## Status: all 6 core scripts + both integration parts covered

TF enrichment, DAR calling, peak fate classification, RNA-seq DEG,
mixed-effects logistic regression, and RNA/ATAC integration Parts 1 and 3
(Part 2 — GSEA-style rank enrichment — is off by the real script's own
default and not covered by either integration kit). This closes out the
full list of scripts identified as central to the manuscript's claims.
