# Test kit: peak fate classification (`01_peak_fate_classification.r`)

Runs your **actual, unmodified** peak-fate script on a small synthetic DAR
dataset with known correct classifications, so a reviewer can verify the
Maintained/Reversed/Converged/Emerged logic (Part A) and the non-vent
Gain/Lose logic (Part B) actually work — without needing the real ATAC-seq
DAR-calling output.

## Why this works without the real pipeline

This script only *requires* two files per contrast: `<stem>_DARs.bed` and
`<stem>_DARs_annotated.txt`, both plain text — no BAMs, no RDS, no
Bioconductor annotation packages needed to satisfy them. Everything else it
can use (called-peak support files, unfiltered `*_FULL_results.txt`,
ChIPseeker/TxDb genome annotation) is optional and degrades gracefully —
confirmed by reading the actual guard conditions in the script (e.g.
`if (exclude_unsupported_emerged && has_called_peaks)`), not assumed.

**The peak-fate script itself is not modified in any way.**

## What's simulated (Part A: BOTv vs BOTCv)

20 synthetic DARs at nc14b, engineered into four known outcomes:

| nc14b DAR indices | Behavior | Expected fate |
|---|---|---|
| 0–7 (8 peaks) | Same coordinates + same-sign log2FC at nc14late | **Maintained** |
| 8–11 (4 peaks) | Same coordinates + opposite-sign log2FC at nc14late | **Reversed** |
| 12–19 (8 peaks) | No corresponding nc14late peak at all | **Converged** |
| — | 5 brand-new nc14late-only peaks (no nc14b overlap) | **Emerged** |

Plus 4 `BOTv_temporal` and 4 `BOTCv_temporal` DARs (some overlapping the
Emerged block) to exercise the `temporal_label` annotation.

## What's simulated (Part B: non-vent, 3 genotype pairs)

Each of the 3 pairs (`BOT`, `BOTC`, `BOTC_oR`) gets 10 temporal DARs (5
gaining / 5 losing accessibility), with 5 of the 10 also overlapping that
genotype's nc14b "anchor" DAR set (`also_nc14b = TRUE`).

## How to install and run

1. Copy this whole `test_kit_peak_fate/` folder to the root of
   `Sun_et_al_2026/` (sibling of `ATACseq_Processing/`) — same layout as the
   TF-enrichment test kit.
2. From inside `test_kit_peak_fate/`, run:
   ```
   Rscript 00_TEST_peak_fate_quickstart.R
   ```

This copies the toy DAR files into `./ATAC_DARS/`, then sources the real
`01_peak_fate_classification.r` script (its own `dar_dir`/`out_dir` config
resolves relative to wherever you run it, so no path editing needed), then
checks the output CSV against the expected counts above.

## What success looks like

```
================ ALL CHECKS PASSED ================
```

with an observed fate table of `Maintained=8, Reversed=4, Converged=8,
Emerged=5` and 30 non-vent rows. Output PDFs land in
`./Overview_Plots/peak_fate/BOTv_vs_BOTCv/` and `.../NonVent/`.

You'll also see `[NOTE]` messages in the console about missing called-peak
files and `*_FULL_results.txt` — **these are expected**, not errors; they're
the script telling you it's skipping the optional rescue/annotation tiers,
exactly as intended for this lightweight test.

## Important: I have not executed this script

Same caveat as the other two test kits — no R/Bioconductor available where
I built this. This one carries more risk than the TF-enrichment kit but
less than the DAR-calling kit: I traced the actual classification logic
(the `case_when()` at line ~1111 for Part A, the direction logic at line
~2597 for Part B) and the graceful-degradation guards by reading the real
code, and cross-checked my expected numbers against that logic by hand —
but I did not trace every line of this 2,771-line script (in particular,
the A6–A14 plotting/statistics sections past the core classification, and
some of the edge-case handling in the reciprocal-overlap join
`best_recip_join()`). **Please run it and paste me the console output**,
especially the `[MISMATCH]` section if anything doesn't match — that'll
tell us exactly where my understanding of the script's logic diverged from
its actual behavior.

## Status: 3 of ~6 pipeline stages covered

TF enrichment, one DAR-calling script, and this peak-fate script. Still
remaining: the other 3 DAR-calling scripts (`02_dar_calling_ventralized_temporal.r`,
`03_dar_calling_nonventralized_temporal_expansion.r`,
`04_mixed_effects_logistic_regression_...r`), the strict-mode peak-fate
variant (`02_peak_fate_classification_strict.r`), the alluvial plotting
scripts, RNA/ATAC integration, and the RNA-seq DEG script. Ask for the next
one when you're ready — and it'd help to know whether the DAR-calling kit
from last time actually ran cleanly, since that result will tell us how
much to trust this same approach going forward.
