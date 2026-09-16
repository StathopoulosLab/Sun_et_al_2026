# Test kit: RNA/ATAC integration (`atac_rna_integration_Step1.R`, Part 1)

Runs your **actual, unmodified** integration script's core concordance
logic — merging ATAC accessibility changes with RNA expression changes by
gene symbol, classifying into concordant/discordant quadrants — on
synthetic data with a known correct merge outcome.

## Scope: Part 1 only, 1 of 5 RNA contrasts

The real script has 3 parts (`--parts` flag, default `"1,3"`). This kit
runs **only Part 1** (`--parts 1`):
- **Part 2** (GSEA-style rank enrichment) is already off by the script's
  own default — not this kit's decision.
- **Part 3** (cluster × fate enrichment) needs `CLUSTER_DIR`, which still
  points at an external path outside this repo (see the `TODO` comment
  in the script's own `CONFIG` section) — not covered by this kit yet.

Of the 5 `RNA_CONTRASTS` the real script defines, this kit supplies toy
data for exactly one: `BOTCv_nc14b_vs_BOTv_nc14b` × its ATAC pair
`BOTCv_vs_BOTv_nc14b` (stem `BOTv_vs_BOTCv`). The other 4 will print
`"Not found"` and get skipped by the script's own `next` in its driver
loop — expected, not a failure.

## Prerequisite fix (already applied, separately)

`DEG_DIR` previously pointed at `../../seRNA-seq_03_19_26/...`, an
external path outside this repo entirely — meaning the real script
couldn't find its RNA input on any machine but the original analyst's.
That's been fixed (changed to match `Analysis_DEG_LimmaBOTvBOTCv.R`'s own
default output path, both scripts now expected to run from the same
working directory). This test kit assumes that fix is already in place.

## What's simulated

60 genes with **identical symbols** in both the toy RNA
`results_full.tsv` (same 8-column format `Analysis_DEG_LimmaBOTvBOTCv.R`
actually produces) and the toy ATAC DAR annotated file — this guarantees
`merge_atac_rna()`'s symbol-based merge succeeds by construction, since
real gene-symbol annotation of ATAC peaks would otherwise require
ChIPseeker + genome annotation packages neither test kit installs.

Engineered as three groups:
- **25 genes**: ATAC and RNA change in the *designed* same direction
  ("concordant" by construction), both clearing their respective
  significance/effect-size thresholds
- **15 genes**: opposite directions ("discordant" by construction), both
  clearing thresholds
- **20 genes**: weak ATAC signal (|log2FC| < 0.50 — won't clear
  `LFC_ATAC`, exercises the "faded, non-significant" point rendering)

## Expected results

```
Merged table: 60 genes
sig_rna:  40/60 (padj < 0.10)
sig_atac: 40/60 (|log2FC| > 0.50)

Quadrant breakdown (exact, computed independently in Python):
  Concordant Open+Up:    19
  Concordant Close+Down: 17
  Discordant Open+Down:  12
  Discordant Close+Up:   12
```

Note the quadrant counts don't map 1:1 onto the "25 concordant / 15
discordant / 20 weak" design groups above — quadrant assignment is by
sign alone (not by significance), so some of the 20 "weak" genes still
land in a concordant or discordant quadrant by chance sign draw. The
**36 total concordant vs 24 discordant** split is the number that
actually matters for validation, and it's independently verifiable from
`test_data`'s generation log.

## How to install and run

1. Copy this whole `test_kit_rna_atac_integration/` folder to the root of
   `Sun_et_al_2026/` (sibling of `ATACseq_Processing/`).
2. From inside `test_kit_rna_atac_integration/`, run:
   ```
   Rscript 00_TEST_rna_atac_integration_quickstart.R
   ```

This copies toy DAR data into `./ATAC_DARS/` and toy DEG output into
`./results/combined_figures/qc_and_deg_botv_botcv/deg_limma/` (matching
the now-fixed `DEG_DIR` default), then calls the real script as a
subprocess with `--parts 1`, then checks the merged output table against
the expected numbers above.

Requires: `GenomicRanges`, `ggplot2`, `dplyr`, `patchwork`, `scales`
(required by the real script); `ggrepel` (optional — gene-label repelling
on the scatter plot, degrades gracefully if absent).

## What success looks like

```
================ PASS ================
60 genes merged, plot + table saved. Core ATAC x RNA concordance logic works.
```

Output lands in `./Overview_Plots/atac_rna_integration/part1_concordance/`
— a combined scatter + quadrant-bar + annotation-feature PDF/PNG, plus the
full merged `.tsv` table.

## Update: confirmed by an actual test run

This kit has been run against the real script and the core merge logic
produced an **exact match** to the predicted numbers: 60/60 genes merged
by symbol, quadrant split `19/17/12/12` (Concordant Open+Up / Concordant
Close+Down / Discordant Open+Down / Discordant Close+Up) — identical to
the independently-computed Python prediction above.

**One real finding along the way, informational rather than a bug**: the
real script writes each quadrant label with a literal embedded newline
(`"Concordant\nOpen+Up"`, used to line-break the legend text in the plot),
and saves the merged table via `write.table(..., quote=FALSE)`. Since
that newline isn't quoted in the output file, it lands as a genuine line
break in the `.tsv`, splitting one data row across two physical lines.
Standard TSV parsers (including R's own `read.table()`) can't cleanly
re-read that file back — confirmed directly when this kit's first version
tried to. This doesn't affect the plots (which is all this script itself
does with the data) but would trip up anyone trying to load
`concordance_*.tsv` into another tool or script later. Worth knowing about
if that file is ever meant to be a reusable output, not just a
plot-generation intermediate.

## Important: I have not executed this script

Same caveat as every other kit — no R available where I built this. I
traced the exact merge logic (`merge_atac_rna()`'s symbol-then-gene_id
matching, the quadrant `case_when()`, the significance flags) and
independently verified the quadrant breakdown in Python, and that's now
confirmed correct by an actual run (see "Update" above). The actual
plotting code (1a/1b/1c) and the driver loop's full behavior across all 5
`RNA_CONTRASTS` were not traced as thoroughly as the merge logic itself —
worth a visual check of the output PDF to confirm the plots themselves
render sensibly.

## Status: 6 of ~6 core pipeline stages covered

TF enrichment, DAR calling, peak fate classification, RNA-seq DEG,
mixed-effects logistic regression, and now RNA/ATAC integration (Part 1).
Remaining, lower-priority items from the original triage: the other 2
DAR-calling script variants, the strict-mode peak-fate script, Part 3 of
this integration script (needs `CLUSTER_DIR` resolved first), and
plotting-only scripts (generally not recommended — see earlier triage).
