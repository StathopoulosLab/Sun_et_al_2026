# Test kit: RNA-seq DEG analysis (`Analysis_DEG_LimmaBOTvBOTCv.R`)

Runs your **actual, unmodified** limma-voom DEG script on a small synthetic
16-embryo counts matrix, exercising every sample-naming convention the
script supports and both of its documented edge-case guards.

## Why this one was easier to verify than the others

This script already takes its input via a `--counts` CLI flag and reads a
plain-text counts matrix — no BAM files, no RDS caches, no Bioconductor
genome-annotation packages. That let me build the toy data directly and
predict the script's sample-classification output (which embryos go to
which group/stage/batch) by exactly replicating its `get_group()`/
`get_stage()`/`EXCLUDE_BCS` regex logic in Python — see the numbers below,
all independently re-derivable from `test_data/counts_matrix.tsv`'s column
names.

**The DEG script itself is not modified in any way.**

## What's in the toy data

16 synthetic embryos across every naming pattern the real script's comments
describe supporting:

| Test case | Sample(s) | What it proves |
|---|---|---|
| Bare/no-stage-tag naming | `bc1` | Defaults correctly to nc14b |
| `foxl1_Dm_*` prefix form | `bc2`, `bc4-6` | Alternate GROUP_A naming matches |
| Bare `BOTv_*` form | `bc3` | Set3 naming convention matches |
| `HLH54F_Dm_*` prefix form | `bc11`, `bc13-16` | Alternate GROUP_B naming matches |
| **B6_BOTC decoy guard** | `bc20`, `bc23` | `foxl1_Dm` + barcode 19-24 correctly excluded (NOT misassigned to FoxL1_BOTv) |
| **EXCLUDE_BCS list** | `bc27` | Explicitly-excluded barcode correctly dropped |
| Batch/sublibrary detection | two `S1`/`S2` prefixes | Batch covariate added, confounding check passes clean |

800 genes, 40 with an engineered genotype effect (log2FC 1.5–3, random
sign) and 40 with an engineered temporal effect (progressing nc14b→nc14d→
gastr), so limma-voom has real signal to detect rather than pure noise.

## Expected results (from `gen_rnaseq_deg_data.py`'s prediction, reproduced here)

```
Total kept embryos: 13 (2 B6_BOTC decoys + 1 EXCLUDE_BCS match correctly dropped from 16)
GROUP_A (FoxL1_BOTv): 6
GROUP_B (HLH54F_BOTCv): 7

Contrast arm sizes:
  1. BOTCv_nc14b_vs_BOTv_nc14b:   A=3 vs B=3
  2. BOTCv_late_vs_BOTv_nc14d:    A=4 vs B=3
  3. BOTCv_late_vs_BOTCv_nc14b:   A=4 vs B=3
  4. BOTv_nc14d_vs_BOTv_nc14b:    A=3 vs B=3
  5. BOTCv_late_vs_BOTv_nc14b:    A=4 vs B=3

All 5 contrasts should run (none SKIP — every arm has >=2 samples).
All 5 should report "2 sublibraries detected — adding batch covariate"
with NO "confounded" warning (batches were deliberately spread across
both conditions in every arm to avoid that edge case).
```

I can't predict the exact DE gene calls (that's real limma/voom statistics
running on the data), but the 40 genotype-effect and 40 temporal-effect
genes should produce non-zero significant hits in the contrasts that
compare across those axes.

## How to install and run

1. Copy this whole `test_kit_rnaseq_deg/` folder to the root of
   `Sun_et_al_2026/` (sibling of `RNAseq_Processing/`).
2. From inside `test_kit_rnaseq_deg/`, run:
   ```
   Rscript 00_TEST_rnaseq_deg_quickstart.R
   ```

This calls the real script as a subprocess with
`--counts ./test_data/counts_matrix.tsv --outdir ./test_output` (its actual
supported CLI flags — see the script's own usage comment at the top), then
checks `summary_table.tsv` against the expected contrast list.

Requires: `limma`, `edgeR` (required by the real script); `ggplot2`,
`ggrepel` (optional — falls back to base-R plotting if absent, per the
script's own `use_ggrepel` check).

## What success looks like

```
[PASS] All 5 expected contrasts completed (none SKIPped).
```

Scroll up in the console output to confirm the sample-classification lines
match the "Expected results" table above. Output lands in
`./test_output/<contrast_name>/` — `results_full.tsv`, `results_sig.tsv`,
and `volcano.pdf`/`.png` per contrast, plus `summary_table.tsv`.

## Important: I have not executed this script

Same caveat as the other three kits — no R available where I built this.
This one carries the **lowest risk of the four so far**: the input format
is plain text I could construct directly, and I independently verified the
sample-classification predictions in Python by replicating the exact regex
logic rather than guessing. The part I couldn't verify by hand is the
actual limma-voom statistics (design matrix construction, `voom()`,
`eBayes()`) — that's real statistical computation I have no way to
hand-check, so if something goes wrong it's more likely to be there than
in the sample-matching logic. Please run it and paste the console output.

## Status: 4 of ~6 pipeline stages covered

TF enrichment, one DAR-calling script, peak fate classification, and this
DEG script. Still remaining: 3 more DAR-calling script variants, the
strict-mode peak-fate script, alluvial plotting, RNA/ATAC integration, and
the RNA-seq preprocessing steps (demux/align/featureCounts — likely similar
in spirit to the ATAC BAM-dependency problem). It'd help to know how the
first three kits actually ran before continuing much further, since that
tells us how reliable this whole approach has been so far.
