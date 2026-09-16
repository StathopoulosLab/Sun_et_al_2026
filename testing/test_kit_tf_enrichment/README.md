# Test kit: TF-enrichment helper functions

A minimal, self-contained test for `ATACseq_Processing/03_TF_Enrichment/00_tf_venn_shared_helpers.R`,
using synthetic peak data instead of the real (large, restricted) sequencing data.

## What's here

```
test_kit/
├── 00_TEST_tf_enrichment_quickstart.R   <- run this
├── test_data/
│   ├── README.md                          <- how the toy data was built + expected numbers
│   ├── DAR_Runt_specific.bed
│   ├── DAR_Cic_specific.bed
│   └── tf_partners/
│       ├── TF_Runt.bed
│       ├── TF_Cic.bed
│       ├── TF_Zld.bed
│       └── TF_Bcd.bed
└── test_output/                           <- created when you run the script
```

## How to install into your repo

Drop this whole `test_kit/` folder in at the root of `Sun_et_al_2026/`
(as a sibling of `ATACseq_Processing/`), so the relative path in the script
(`../ATACseq_Processing/03_TF_Enrichment/00_tf_venn_shared_helpers.R`)
resolves correctly:

```
Sun_et_al_2026/
├── ATACseq_Processing/
├── RNAseq_Processing/
├── Runt_ChIPseq/
└── test_kit/          <- new
```

## How to run it

```
cd Sun_et_al_2026/test_kit
Rscript 00_TEST_tf_enrichment_quickstart.R
```

Requires the same R packages `00_tf_venn_shared_helpers.R` already needs:
`GenomicRanges`, `ggplot2`, `dplyr`, `tidyr`, `stringr`, `scales`,
`VennDiagram`, `futile.logger`, `patchwork`.

## What success looks like

The script ends by printing:
```
================ ALL TESTS PASSED ================
```

and writes three output files to `test_output/` that you can open to
visually confirm the plots render correctly:
- `Explanatory_Runt_specific_panel.pdf`
- `Explanatory_Cic_specific_panel.pdf`
- `TEST_venn_runt_vs_cic.pdf`

If any `stopifnot()` check fails, the script stops immediately and tells you
which expected value didn't match — see `test_data/README.md` for exactly
how each expected number was derived, so you can trace a mismatch back to
whichever function produced it.

## Important: I have not executed this script

I built and hand-verified the toy data and expected numbers in Python
(no R available in the environment I was working in), and wrote this script
by reading `00_tf_venn_shared_helpers.R`'s actual function signatures —
but I could not install the Bioconductor/R dependencies to run it myself.
**Please run it once locally before relying on it** — if `library(GenomicRanges)`
or any other `library()` call fails, that just means a package needs
installing, not that the test itself is wrong. If you hit any error beyond
a missing package, paste it to me and I'll help debug it.

## Status: 1 of ~6 pipeline stages covered

This covers `03_TF_Enrichment/00_tf_venn_shared_helpers.R` only — the most
self-contained, function-based part of the pipeline (no BAM files or
hardcoded sample lists to worry about). The DAR-calling scripts
(`02_DAR_Calling/`), peak-fate classification (`05_Peak_Fate/`), and RNA-seq
DEG analysis will need a similar treatment, but each requires first
reading through hardcoded file paths/sample sheets to figure out what
minimal parameterization (if any) is needed to make them toy-data-runnable.
Ask for the next stage when you're ready.
