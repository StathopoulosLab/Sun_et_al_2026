# RNA-seq analysis pipeline (single-embryo CEL-Seq2 / seRNA-seq)

Processes pooled, barcoded single-*Drosophila*-embryo RNA-seq (CEL-Seq2-style)
from raw FASTQs through demultiplexing, alignment, UMI counting, QC filtering,
normalization, batch correction, and downstream differential expression /
classification / visualization.

All scripts assume they are run from the repository root, with a shared
`config.py` in the same directory defining:

```python
OUTPUT_ROOT      = "results"          # root output directory
STAR_INDEX       = "/path/to/STAR_index"
GTF_FILE         = "/path/to/dm6.ensGene.gtf"
GENE_LENGTHS     = "dm6_gene_lengths.tsv"
THREADS_DESKTOP  = 4
THREADS_HPC      = 16
```

(This file is local to each install and is not included here — create your
own before running any step.)

## Pipeline order

| Step | Script | Purpose | Key output |
|---|---|---|---|
| 0 | `Step00_RNAseq_GeneLengths.py` | One-time: build a gene-length table from the GTF | `dm6_gene_lengths.tsv` |
| 1 | `Step01_RNAseq_Demux.py` | Demultiplex a pooled run into per-embryo FASTQs by barcode | `results/<sample>/demux/*_R1.fastq`, `*_R2.fastq` |
| 2 | `Step02_RNAseq_Align.py` | Align each embryo's R2 reads with STAR | `results/<sample>/mapped/*.sorted.bam` |
| 3 | `Step03_RNAseq_FeatureCounts.py` | Tag reads with their gene (featureCounts `-R BAM`) | `results/<sample>/mapped/*.sorted.featureCounts.bam` |
| 4 | `Step04_RNAseq_UMICollapse.py` | Collapse UMIs per gene per embryo | `results/<sample>/counts/*_counts.txt` |
| 5 | `Step05_RNAseq_MergeCounts.py` | Merge one sample's per-embryo counts into a matrix | `results/<sample>/counts/<sample>_counts_matrix.tsv` |
| 5B | `Step05B_RNAseq_MergeSamples.py` | Merge all samples/sets into one genes × embryos matrix, with group/stage/set metadata | `results/combined/<prefix>_counts_matrix.tsv`, `..._sample_groups.tsv` |
| 5C | `Step05C_RNAseq_QCFilter.py` | Group- and stage-aware QC filtering (UMI/gene thresholds, mito%, diagnostic-gene zero-expression check) | `results/combined/<prefix>_counts_matrix_qc.tsv` |
| 6 | `Step06_RNAseq_TPMNormalization.py` | UMI counts → TPM, optional quantile normalization | `results/combined/<prefix>_tpm_matrix.tsv` (+ QN variants) |

`Step05B_..._legacy.py` and `Step05C_..._legacy.py` are earlier versions of
those two steps, kept for provenance/reproducibility of intermediate results
generated before the current versions were adopted. Use the non-`_legacy`
version for new runs.

### Example end-to-end run for one sample

```bash
python Step01_RNAseq_Demux.py --r1 raw/sample_4embryos_R1.fastq.gz \
                               --r2 raw/sample_4embryos_R2.fastq.gz \
                               --barcodes barcodes/4embryos.tsv

python Step02_RNAseq_Align.py --sample sample_4embryos --hpc
python Step03_RNAseq_FeatureCounts.py --sample sample_4embryos --hpc
python Step04_RNAseq_UMICollapse.py --sample sample_4embryos
python Step05_RNAseq_MergeCounts.py --sample sample_4embryos
```

### Combining samples and normalizing

```bash
python Step05B_RNAseq_MergeSamples.py \
    --wt-sample WT_seRNAseq_4embryos_set1 \
    --mutant-sample Mutant_seRNAseq_19embryos_set2 \
    --set2-sample 28232_28233_28271_28272_28273_S1 \
    --prefix all_genotypes

python Step05C_RNAseq_QCFilter.py --combined all_genotypes --min-umis 8000 --min-genes 5000 --goi-check
python Step06_RNAseq_TPMNormalization.py --combined all_genotypes --qc --quantile
```

## Support / config modules (imported, not run directly)

| File | Purpose |
|---|---|
| `Config_SampleMetadata.py` | Single source of truth for genotype vocabulary, barcode→genotype maps, validated embryo lists, manual exclusions, and QC thresholds. Imported as `Config_SampleMetadata` (commonly aliased `SC`) by most downstream scripts. |
| `Config_DataLoader.py` | Standard loader (`DL.load()`) that returns RUVg-corrected, QC-filtered, quantile-normalized data with group/stage labels attached, used by all plotting/classification scripts. |

Embryos are dropped from analysis for two independent reasons, both tracked
in `Config_SampleMetadata.py`:
- `MANUAL_EXCLUDE` — embryos failing basic QC (demux failure, sub-threshold
  UMI capture, ambiguous barcode). Applied everywhere via
  `apply_manual_exclude()` / `Step05C_RNAseq_QCFilter.py`.
- Per-script `EXCLUDE_EMBRYOS` (in the plotting/classification scripts) —
  a small, secondary set of outliers excluded only from specific
  visualizations, not part of the primary QC pipeline.

## Downstream analysis / plotting (run after Step06)

| File | Purpose |
|---|---|
| `Analysis_RUVgNormalization.R` | RUVg batch-effect correction on the combined count matrix. |
| `Analysis_DEG_LimmaBOTvBOTCv.R` | limma differential expression between BOTv and BOTCv genotypes. |
| `Plot_PCA_UMAP.py` | PCA / UMAP embedding and QC-labeled variants. |
| `Plot_Heatmaps_Informative.py` | Heatmaps of informative/marker genes by genotype and call. |
| `Plot_Heatmaps_DevelopmentalCascade.py` | Heatmaps spanning all developmental timepoints. |
| `Plot_PositionalClassification.py` | Classifies embryos by anterior-posterior / dorsal-ventral marker gene expression; produces violin plots, heatmaps, and per-embryo classification tables. |

These all import `Config_SampleMetadata` and `Config_DataLoader`, so both
must be present in the same directory (or on `PYTHONPATH`).

## Dependencies

- **Python**: `pandas`, `numpy`, `pysam` (Step04 only)
- **R**: `limma`, `edgeR`/`RUVSeq` (for `Analysis_RUVgNormalization.R`), plus
  standard Bioconductor annotation packages for dm6
- **External tools**: STAR, featureCounts (Subread), and (for
  `Step00_RNAseq_GeneLengths.py`) a dm6 GTF annotation

## Notes

- All scripts read paths through `config.py` and `Config_SampleMetadata.py`
  rather than hardcoding them — set `OUTPUT_ROOT` etc. once in `config.py`
  and everything downstream resolves consistently.
- A one-off patch that updated `Config_DataLoader.py`'s embryo-group
  assignment logic (dropping unvalidated "unknown" embryos rather than
  guessing their group) has already been applied to the version of
  `Config_DataLoader.py` in this repo; there is no separate patch script.
