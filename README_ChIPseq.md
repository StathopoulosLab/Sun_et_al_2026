---
editor_options: 
  markdown: 
    wrap: 72
---

# ChIP-seq preprocessing

Paired-end ChIP-seq preprocessing for *Drosophila* embryos, from raw
FASTQs through alignment, filtering, peak calling, and normalized
signal-track (bigWig) generation, run as SLURM (`.sub`) jobs on an HPC
cluster.

## Setup

Every script has a small config block near the top:

``` bash
PROJECT_DIR="/path/to/your/project"
CONDA_BASE="$(conda info --base)"
```

Set `PROJECT_DIR` to your own project root (expects
`ChIP_seq/ChIP_data/` and `reference_genome/` subdirectories beneath
it). All scripts activate a conda environment named `atac_env` (shared
with the ATAC pipeline) — it needs `bowtie2`, `samtools`, `macs2`,
`deeptools`, `bedtools`, and UCSC tools (`bedGraphToBigWig`) installed.

Sample and control BAM paths (e.g. `Dm_ChIP_D_GFPAb_i6_JS`,
`Dm_ChIP_D_INPUT_i19_JS`, `Dm_ChIP_D_IgG_i12_JS`) are hardcoded inside
`Step02*` and `Step03*` scripts rather than passed as arguments — edit
the sample-definition block near the top of each script for your own
samples.

## Pipeline

```         
Step01_ChIPseq_IR_Trimmming-Mapping_2_.sub
        │  Adapter/quality trimming + Bowtie2 alignment to dm6
        ▼
Step01B_ChIPseq_Enhanced_Filtering_1_.sub          (per-sample, incl. controls)
        │  Additional filtering for cleaner tracks (dedup, blacklist, quality)
        ▼
Step01C_Merge_Filtered_Controls.sub                (controls only)
        │  Merge enhanced-filtered INPUT/IgG replicates into pooled controls
        │
        ├── BASELINE peak-calling branch ─────────────┐
        │   (uses standard dedup BAMs from Step01)     │
        │                                              ▼
        │                  Step02_ChIPseq_IR_PeakCalling_..._igg-input.sub
        │                  MACS2, multiple stringencies, vs INPUT and/or IgG
        │
        └── FILTERED peak-calling branch ──────────────┐
            (uses enhanced-filtered BAMs from           │
             Step01B / Step01C)                         ▼
                              Step02B_ChIPseq_PeakCalling_Filtered_1_.sub
                              Same MACS2 logic as Step02, cleaner input BAMs

Signal-track (bigWig) generation — three alternative methods, all
consuming the filtered BAMs (Step01B/Step01C) or standard BAMs (Step01):

  Step03B_ChIPseq_Enhanced_Tracks_ATAC_style.sub
      ATAC-style background normalization, single control

  Step03B_ChIPseq_Enhanced_Tracks_ATAC_style_wMerg_Igg_1_.sub
      Same as above, using the merged IgG control from Step01C

  Step03C_ChIPseq_MACS2_bdgcmp_Tracks_1_.sub
      MACS2 bdgcmp-based tracks ("paper method")
```

Peak calling (Step02/Step02B) and track generation (Step03B/Step03C) are
independent downstream uses of the filtered BAMs — run whichever
branch(es) you need; they don't depend on each other.

## Script details

| Script | Role | Input | Output |
|------------------|------------------|------------------|------------------|
| `Step01_ChIPseq_IR_Trimmming-Mapping_2_.sub` | Step 1 | raw paired FASTQs | trimmed FASTQs, aligned/dedup BAM |
| `Step01B_ChIPseq_Enhanced_Filtering_1_.sub` | Step 1B | dedup BAM from Step01 | enhanced-filtered BAM |
| `Step01C_Merge_Filtered_Controls.sub` | Step 1C | filtered INPUT/IgG replicate BAMs | pooled control BAM(s) |
| `Step02_ChIPseq_IR_PeakCalling_1_14_26_specific_V2_allstringency_igg-input.sub` | Step 2 (baseline) | standard dedup ChIP + INPUT/IgG BAMs | MACS2 peaks, multiple stringencies |
| `Step02B_ChIPseq_PeakCalling_Filtered_1_.sub` | Step 2B (filtered) | enhanced-filtered ChIP + pooled control BAMs | MACS2 peaks, multiple stringencies |
| `Step03B_ChIPseq_Enhanced_Tracks_ATAC_style.sub` | Step 3B (track method 1) | filtered ChIP + control BAM | bigWig (ATAC-style normalized) |
| `Step03B_ChIPseq_Enhanced_Tracks_ATAC_style_wMerg_Igg_1_.sub` | Step 3B variant | filtered ChIP + merged IgG BAM | bigWig (ATAC-style normalized) |
| `Step03C_ChIPseq_MACS2_bdgcmp_Tracks_1_.sub` | Step 3C (track method 2) | filtered ChIP + control BAM | bigWig (MACS2 bdgcmp) |

## Example run for one sample

``` bash
sbatch Step01_ChIPseq_IR_Trimmming-Mapping_2_.sub BASENAME
sbatch Step01B_ChIPseq_Enhanced_Filtering_1_.sub BASENAME
# run Step01B on each control replicate too, then:
sbatch Step01C_Merge_Filtered_Controls.sub

# Peak calling — edit the sample block at the top of one of these, then:
sbatch Step02_ChIPseq_IR_PeakCalling_1_14_26_specific_V2_allstringency_igg-input.sub
# or, using filtered BAMs:
sbatch Step02B_ChIPseq_PeakCalling_Filtered_1_.sub

# Tracks — pick one method, edit its sample block, then:
sbatch Step03B_ChIPseq_Enhanced_Tracks_ATAC_style.sub
sbatch Step03C_ChIPseq_MACS2_bdgcmp_Tracks_1_.sub
```
