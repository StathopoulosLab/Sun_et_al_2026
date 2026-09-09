# ATAC-seq preprocessing (fastq → bam → bigwig)

Paired-end ATAC-seq preprocessing for *Drosophila* embryos, from raw FASTQs
through alignment, filtering, peak calling, and normalized signal-track
(bigWig) generation, run as SLURM (`.sub`) jobs on an HPC cluster.

## Setup

Every script has a small config block near the top:

```bash
PROJECT_DIR="/path/to/your/project"
CONDA_BASE="$(conda info --base)"
```

Set `PROJECT_DIR` to your own project root (expects `ATAC_seq/ATAC_data/`,
`ATAC_seq/scripts/`, and `reference_genome/` subdirectories beneath it) and
either let `CONDA_BASE` auto-resolve or hardcode it. All scripts activate a
conda environment named `atac_env` — create one with `bowtie2`, `samtools`,
`macs2`, `deeptools` (for `bamCoverage`), `bedtools`, and UCSC tools
(`bedGraphToBigWig`) installed.

## Pipeline

```
Step01_ATAC_FastQC_TrimGalore_IR.sub
        │  FastQC + TrimGalore (adapter/quality trimming)
        ▼
Step02_ATAC_Align2dm6-filt_IR_fixed.sub
        │  Bowtie2 alignment to dm6 (--very-sensitive, -X 700)
        ▼
Step03_ATAC_Post-Align2dm6-filt_IR_2_.sub
        │  samtools fixmate → markdup (dedup) → chrM removal
        │
        ├──────────────────────────────┬─────────────────────────────┐
        ▼                               ▼
Step04_ATAC_MACS2_Narrow-Broad_        Step04A_ATAC_Norm_IR_merged_3_.sub
Peakcalling.sub                        │  bamCoverage, CPM normalization
   (peak calling branch)               ▼
   narrow + broad peaks               Step04B_ATAC_Norm2_IR_merged_3_.sub
                                       │  background/z-score normalization,
                                       │  blacklist removal, bedGraphToBigWig
                                       ▼
                                   final .bigWig
```

After Step03, the pipeline branches into two independent uses of the same
filtered BAMs:

- **Peaks**: `Step04_ATAC_MACS2_Narrow-Broad_Peakcalling.sub` calls narrow
  and broad peaks directly with MACS2.
- **Signal tracks**: `Step04A` → `Step04B` produce normalized, blacklist-
  filtered bigWig tracks for genome browser visualization. These two steps
  must be run in order — Step04B consumes Step04A's output.

## Script details

| Script | Input | Output |
|---|---|---|
| `Step01_ATAC_FastQC_TrimGalore_IR.sub` | raw paired FASTQs (`<BASENAME>_R1/R2.fastq.gz`) | FastQC reports, trimmed FASTQs |
| `Step02_ATAC_Align2dm6-filt_IR_fixed.sub` | trimmed FASTQs | `<BASENAME>_aligned.sam` |
| `Step03_ATAC_Post-Align2dm6-filt_IR_2_.sub` | Tn5-shifted, blacklist-filtered BAM | deduped, chrM-removed BAM (`*_noq_rmdup.noChrM.bam`) + mapping stats |
| `Step04_ATAC_MACS2_Narrow-Broad_Peakcalling.sub` | filtered BAMs from Step03 | MACS2 narrow/broad peak sets |
| `Step04A_ATAC_Norm_IR_merged_3_.sub` | filtered BAMs from Step03 | CPM-normalized bedgraph (`bamCoverage`) |
| `Step04B_ATAC_Norm2_IR_merged_3_.sub` | Step04A bedgraph | blacklist-filtered, background-normalized `.bigWig` |

## Example run for one sample

```bash
sbatch Step01_ATAC_FastQC_TrimGalore_IR.sub sample_set BASENAME
sbatch Step02_ATAC_Align2dm6-filt_IR_fixed.sub BASENAME
sbatch Step03_ATAC_Post-Align2dm6-filt_IR_2_.sub BASENAME

# Peaks branch
sbatch Step04_ATAC_MACS2_Narrow-Broad_Peakcalling.sub BASENAME

# Bigwig branch (run in order)
sbatch Step04A_ATAC_Norm_IR_merged_3_.sub BASENAME
sbatch Step04B_ATAC_Norm2_IR_merged_3_.sub BASENAME
```
