# Test kit: mixed-effects logistic regression (`04_mixed_effects_logistic_regression_state_selective_accessibility.r`)

Runs your **actual, unmodified** mixed-effects script — real `lme4`/`lmerTest`
model fitting, binary GLMM, optional `glmmLasso` regularization, and
bootstrap stability checks — on synthetic data, for two of its ~27
contrasts.

## Update: confirmed by an actual test run

**The `"not using"` path was NOT a live blocker** — all 19 TF datasets
loaded successfully in an actual test run. Still worth a quick sanity
check on your real data (the folder name is a strange thing to have in a
production path), but it's confirmed not to crash the loading step.

**A real, different problem showed up in the actual test run**: every
single model fit failed, identically, across both contrasts and all 4
strata (8 for 8). The real script swallows the actual `lme4` error message
(`tryCatch(..., error = function(e) NULL)` at line ~1320 — worth flagging
to whoever maintains this: silently discarding the error text makes real
debugging much harder than it needs to be), so the console only showed
`"Fit failed."` with no detail.

**Diagnosis**: this was a toy-data design bug, not a real-script bug. The
model formula is `logFC ~ TF1 + ... + TF19 + (1|chr)` — a random-effects
term needs more than one level to estimate variance from, and the original
toy data put every peak on a single chromosome (`chr2R`). `lme4` cannot
fit `(1|chr)` when `chr` never varies. **Fixed**: toy peaks are now spread
across 5 chromosomes (`chr2L`/`chr2R`/`chr3L`/`chr3R`/`chrX`, 24 each), and
scaled up from 40 to 120 total peaks per contrast (a healthier
observations-to-TF-predictors ratio — 19 TF fixed effects fit against only
~30 peaks at the strictest threshold was thin even setting the
single-chromosome issue aside). This fix has not yet been re-run — please
try again with this updated version.

## Prerequisite fix (already applied, separately)

This script previously hard-`stop()`ped on a missing
`DAR_ULTIMATE_COMPLETE_LOADER.r` file that isn't in the repo. That's been
fixed already (the dead `source()` calls were removed — confirmed the
script never actually used anything from that framework). This test kit
assumes that fix is already in place.

## Scope: 2 of ~27 contrasts

The real script defines contrasts across two genotype splits
(`split1_configs`, `split2_configs`). This kit supplies toy DAR data for
exactly two, chosen to exercise both splits and the stem-specific
sign-negation logic:

| Contrast | Split | In `NEGATE_LFC_STEMS`? |
|---|---|---|
| `BOT_vs_BOTR` | 1 (non-ventralized) | No |
| `BOTv_vs_BOTCv` | 2 (ventralized) | **Yes** — script negates this stem's logFC sign internally |

The other ~25 contrasts will print `[missing]` for their DAR file and get
skipped by the real script's own `if (length(gr) == 0) SKIP` check in
`run_and_save_lmm()` — **this is expected**, not a failure. Building toy
data for all 27 wasn't practical; these two exercise the shared core logic
(TF design matrix construction, threshold stratification, model fitting)
that every other contrast also goes through.

## What's simulated

**DAR data**: both contrasts reuse the same 120 genomic loci, spread
across 5 chromosomes (`chr2L`/`chr2R`/`chr3L`/`chr3R`/`chrX`, 24 each —
necessary so the model's `(1|chr)` random-effect term has real variance to
estimate; see "Update" above), 200bp peaks at 1500bp spacing within each
chromosome, with independently assigned `log2FoldChange` per stem:
- 45 peaks: strong positive (0.6 to 1.5)
- 45 peaks: strong negative (-1.5 to -0.6)
- 30 peaks: weak (±0.1 to 0.2)

This guarantees **90 peaks clear every |logFC| threshold** (0.50, 0.35,
0.25) and all 120 are available for the unfiltered `all_peaks` tier —
comfortably clearing `run_lmm_stratum()`'s own "fewer than 20 peaks" skip
floor at every stratum.

**TF partner peaks**: all 19 files at the real script's exact hardcoded
paths (these are **required**, not optional — `import()` has no
missing-file guard). Each TF binds a random 25–95 of the 120 loci, giving
every TF column real variance so none get dropped by the "no variable TF
columns" check.

## How to install and run

1. Copy this whole `test_kit_mixed_effects/` folder to the root of
   `Sun_et_al_2026/` (sibling of `ATACseq_Processing/`).
2. From inside `test_kit_mixed_effects/`, run:
   ```
   Rscript 00_TEST_mixed_effects_quickstart.R
   ```

This copies the toy TF/DAR data into place, sources the real script, then
checks that both contrasts produced `lmm_coefficients.txt` +
`lmm_result.rds` for all 4 strata (3 thresholds + `all_peaks`).

Requires: `lme4`, `lmerTest`, `GenomicRanges`, `rtracklayer`, `ggplot2`,
`dplyr`, `tidyr` (all required by the real script). `glmmLasso` is
optional — if installed, its regularized models run too; if not, the real
script prints a warning and skips them gracefully.

## What success looks like

```
[PASS] Both toy contrasts produced coefficient + RDS output for all 4 strata.
```

Expect a LOT of console noise — `[missing]` lines for the ~25 untested
contrasts, `SKIP (0 DARs)` lines, and (if `RUN_BOOTSTRAP_STABILITY` is on,
which it is by default) `[bootstrap] running 200 resamples...` for each of
the 8 strata (2 contrasts × 4 strata) that do have data. That's all
expected. Output lands in dated folders (`Split1_NonVent_newcic_<date>/`,
`Split2_Vent_newcic_<date>/`) since the real script names its output
directories with today's date by design.

## Important: I have not executed this script

Same caveat as every other kit — no R available in the environment I built
this in. This is the **most complex script tested so far** — real mixed
models (`lme4::lmer`/`lmerTest`), a collinearity-flagging step, a
bootstrap resampling loop, and optional `glmmLasso` regularization with
BIC-guided lambda selection. I verified the toy data clears every
documented size/variance floor in the code (n≥20 peaks per stratum, TF
column variance > 0), and confirmed the two contrast keys/stems exactly
match the real script's `split1_configs`/`split2_configs` definitions —
but I have not traced the full mixed-model fitting internals
(`run_lmm_stratum`'s remaining ~600 lines, `fit_glmmlasso_safe`,
`bootstrap_stability_check`) line-by-line the way I did for the simpler
scripts. **Please run this and paste the full console output** — given the
complexity, I'd genuinely expect a real chance of something breaking here
that I didn't anticipate.

## Status: 5 of ~6 core pipeline stages covered

TF enrichment, one DAR-calling script, peak fate classification, RNA-seq
DEG, and now mixed-effects logistic regression. RNA/ATAC integration
(`atac_rna_integration_Step1-3.R`) is the remaining high-priority one —
likely underlies whichever figure combines chromatin and expression
findings.
