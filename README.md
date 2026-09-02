# Sun_et_al_2026

# Temporal remodeling and regulatory plasticity resolve adjacent mesodermal states

Overview

This repository contains the preprocessing and analysis pipelines used to generate the ATAC-seq, ChIP-seq, and single-embryo RNA-seq results in the paper. We profiled chromatin accessibility and transcription factor occupancy across a panel of Drosophila melanogaster genotypes spanning wild type, BOT-background (bcd⁻ osk⁻ tsl⁻), and BOTC-background (BOT + cic⁻) embryos, with and without toll10b-driven ventralization, at nuclear cycle 14 (early and late) and gastrulation stages, alongside matched single-embryo RNA-seq to quantify transcriptional output.


# Embryo segmentation and expression analysis (MATLAB)
embryo_segmentation_and_expression_analysis_ESEA_V1.m segments embryos from raw microscopy images and extracts intensity profiles for expression analysis.

embryo_segmentation_and_expression_analysis_ESEA_V1.m
Version requirements:
MATLAB: 2021
Data Format:
tif; lif; czi; mat
Briefly: This code accepts an image file, converts to uint16 and then accepts a user defined threshold to segment the embryo.
The user then picks the embryo mask and the intensity profiles are generated.


Contact

