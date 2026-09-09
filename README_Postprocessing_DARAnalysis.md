# Drosophila ATAC-seq / RNA-seq Chromatin Fate-Decision Pipeline

Analysis code for a multi-genotype ATAC-seq (+ RNA-seq integration) study of
chromatin accessibility dynamics across a *Drosophila* embryo genotype series,
built around a Cic-gated, Runt-dosage-tunable fate-decision model.

This repo holds the analysis scripts only — no sequencing data, BAMs,
BigWigs, or count matrices are included. Every script expects a local
project directory laid out roughly as described in [Data layout](#data-layout)
below, and reads/writes relative paths from wherever it's run.

## Pipeline stages

Scripts are grouped by stage; the numeric prefixes reflect run order, not
strict dependency (some stages, e.g. QC, can run anytime once BAMs/peaks
exist).

### `01_peak_universe/`
| Script | Purpose |
|---|---|
| `01_merge_replicate_peaks_to_group_beds.sh` | Per-genotype-group union and reproducible (≥2 replicate) merged peak BEDs from MACS2 narrowPeak calls. |
| `02_build_genotype_peak_universe.sh` | Builds the full/ventralized/non-ventralized peak universes used as the counting substrate for DAR calling. |

### `02_dar_calling/`
| Script | Purpose |
|---|---|
| `01_dar_calling_seven_genotype_nc14b.r` | Core seven-genotype nc14b DAR analysis (DESeq2 + limma-voom + RUVg), all pairwise/background contrasts. |
| `02_dar_calling_ventralized_temporal.r` | DAR calling for the ventralized (BOTv/BOTCv) genotype pair across nc14b → nc14late, plus derived temporal categories. |
| `03_dar_calling_nonventralized_temporal_expansion.r` | Temporal expansion of the non-ventralized genotype series to later timepoints. |
| `04_mixed_effects_logistic_regression_state_selective_accessibility.r` | Mixed-effects logistic regression modeling of state-selective accessibility: `lme4::glmer` (binomial family, logit link, `(1\|chr)` random intercept) per stage, classifying DARs by BOTv- vs BOTCv-preferential accessibility from TF-occupancy features, with a bootstrap stability check. Significance is from asymptotic Wald *z*-tests on a single full-data fit — see [Note on lmerTest](#note-on-lmertest) below. **This is the version used in the paper.** |

These three DAR-calling scripts are genuinely different pipelines, not
copies of one another with a label swapped — the genotype panels, replicate
counts, and RUVg control-selection logic differ by design between the
ventralized and non-ventralized series. They are not merged into a single
parameterized script for that reason; see [Note on the vent /
non-vent scripts](#note-on-the-vent--non-vent-scripts) below.

### `03_tf_enrichment/`
| Script | Purpose |
|---|---|
| `00_tf_venn_shared_helpers.R` | Shared config + helpers sourced by the Venn/heatmap scripts: master genotype/concept color palettes, `make_venn2()`/`make_venn3()`, the 17-TF partner-peak panel, and the TF "explanatory power" engine. Already condition-agnostic — used by both splits. |
| `01_venn_overlap_and_tf_enrichment.r` | Venn overlap + TF-enrichment analysis for one condition's DAR sets. `CONDITION_NAME`/`CONDITION_TAG` at the top select which condition; ships configured for Split 2 (ventralized). |
| `02_tf_stratified_accessibility_heatmaps.R` | deepTools heatmaps/profiles, TF-co-binding-stratified (Runt-bound / Cic-bound / both / unbound). `FLANK_BP` at the top sets the computeMatrix flanking window — defaults to **2000 (±2kb)**; set to 5000 for the ±5kb variant used elsewhere in the paper. |

### `04_plotting/`
| Script | Purpose |
|---|---|
| `01_peak_fate_classification.r` | Full peak-fate classification (Maintained / Reversed / Converged / Emerged) — Part A: BOTv vs BOTCv nc14b→nc14late; Part B: per-genotype non-ventralized temporal fate. |
| `02_peak_fate_classification_strict.r` | Strict-mode re-derivation of Part A's fate table (no trend/imputation/union/depth/HCR rescue). Used for the paper's strict-criteria sensitivity check. |
| `03_peak_fate_alluvial_plots.r` | Alluvial diagrams of peak fate for the default (non-strict) classification. |
| `04_peak_fate_alluvial_plots_strict.r` | Same alluvial plots, reading from `02_peak_fate_classification_strict.r`'s output instead. |
| `05_genome_browser_coverage_tracks.py` | Genome-browser-style coverage track plotting from BigWigs, with a `--mode` flag (`all`, `ventralized`, `non_ventralized`, `*_timepoints`, `CVM`, `CVM_solo`) selecting which genotype/timepoint panel to draw. |

`05_genome_browser_coverage_tracks.py` is the template for how the other
plotting scripts were generalized: one script, a runtime flag/config block
selects the genotype panel, rather than a separate near-duplicate file per
condition.

### `05_rna_integration/`
| Script | Purpose |
|---|---|
| `atac_rna_integration_v1.R` | v1 — original ATAC/RNA integration pipeline. |
| `atac_rna_integration_v2.R` | v2 — revised CRM/DAR matrix loading. |
| `atac_rna_integration_v3.R` | v3 — current version: label-clipping fixes, refined thresholds. |

All three versions are included (not just the latest) since different
figures/analyses in the paper draw on different versions' outputs.

### `06_qc/`
| Script | Purpose |
|---|---|
| `01_frip_qc.sh` | Local FRiP (fraction of reads in peaks) calculator against already-called MACS2 peaks and local BAMs. |

## Data layout

Scripts assume (and will create as needed) a project directory with
subfolders along these lines — none of this is included in the repo:

```
project_root/
├── ATAC_NarrowPeaks/              # MACS2 narrowPeak calls, one per sample
├── Merged_Peak_BEDs/              # output of 01_merge_replicate_peaks_to_group_beds.sh
├── Control_Bams/                  # BAMs, *_noq_rmdup.noChrM.bam
├── counts_cache/                  # cached featureCounts matrices
├── Generate_fresh_counts/Output/  # per-contrast DAR output (02_dar_calling/*)
├── BigWig_Files_CPM/              # per-sample CPM bigWigs
├── BigWig_Diagnostic_All/merged/  # merged/background-subtracted bigWigs
├── data/chipseq/                  # Runt ChIP-seq BED
├── data/tf_partners/              # Cic + partner-TF ChIP-seq BEDs
└── Overview_Plots/                # figure output root
```

Each script also expects an optional project-local path config
(`paths_config.r` / `NV_PATHS.r` in the original project) that isn't part of
this repo — it just held machine-specific absolute paths. Point the
`*_dir`/`*_base` variables near the top of each script at your own layout,
or write a small `paths_config.r` that sets them before sourcing.

## Note on the vent / non-vent scripts

Several scripts were originally split into near-duplicate "ventralized" and
"non-ventralized" copies that differed only in which genotype panel, output
directory, and title string they used — not in the underlying logic. Where
that was true (the Venn/TF-enrichment script and the heatmap script), the
two copies were consolidated into one script with a small `CONDITION_NAME`
/ `CONDITION_TAG` (or `FLANK_BP`, for the heatmap window size) configuration
block at the top, so the same script serves either condition depending on
what you point it at.

Where the "vent" and "non-vent" scripts actually encode different biology —
different genotype panels, different replicate structures, different tier
logic for the DAR-calling engines in `02_dar_calling/` — they were kept as
separate scripts rather than force-merged. Collapsing those into a single
parameterized script would risk silently changing the statistics behind
published results, so only cosmetic/labeling cleanup was applied there
(hardcoded hardware/local-machine comments removed; everything else
untouched).

## Note on lmerTest

`04_mixed_effects_logistic_regression_state_selective_accessibility.r`
loads both `lme4` and `lmerTest`, but only the continuous-model code path
(not present in the paper's reported analysis) calls `lmerTest::lmer()`
directly. The binary/state-selective model reported in the paper calls
`glmer()` from base `lme4`, and `lmerTest` does not override `summary()`
for `glmer` fits — only for `lmer` fits. So the significance calls actually
used (coefficients, model-based SE, and two-sided *P* values) are the
standard asymptotic Wald *z*-tests from `lme4::glmer`, not a `lmerTest`
Satterthwaite/KR-df result. Methods text and any Key Resources Table entry
describing the binary-model significance calls should cite `lme4::glmer`
only.

## Requirements

R (≥4.x) with: `DESeq2`, `limma`, `edgeR`, `RUVSeq`, `GenomicRanges`,
`rtracklayer`, `ggplot2`, `dplyr`, `tidyr`, `patchwork`, `scales`,
`ggalluvial`, `VennDiagram`, `stringr`, `lme4`/`lmerTest`, `glmmLasso`.

Python (≥3.8) with: `pyBigWig`, `matplotlib`, `numpy`, `pandas`.

Command-line: `bedtools`, `samtools`, `deepTools` (`computeMatrix`,
`plotHeatmap`, `plotProfile`).

## Provenance

These scripts were sanitized for public release from an internal lab
repository: local absolute paths, machine-specific comments, and
internal-only notes were removed or genericized. Scientific content
(genotype definitions, thresholds, statistical methods, bug-fix notes) was
left intact as methods documentation.
