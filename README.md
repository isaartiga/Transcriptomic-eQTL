# Transcriptomic-eQTL
## Overview
This repository contains the scripts, examples of data, and results for an eQTL analysis conducted as part of a Master's Thesis. The study focuses on investigating the regulatory effects of genetic variants on gene expression in a cohort of Relapsing-Remitting Multiple Sclerosis (RRMS) patients and healthy controls.

## Objectives
- To explore transcriptional differences between RRMS patients and healthy controls by identifying differentially expressed genes (DEGs).
- To perform stringent quality control (QC) procedures on genotyping data to ensure the accuracy and reliability of genetic variants and sample integrity for cis-eQTL analyses.
- To conduct cis-eQTL analysis using DEGs to uncover genetic variants influencing gene expression in pDCs.
- To compare the genotypic variability of key cis-eQTL-associated single nucleotide polymorphisms (SNPs) between RRMS patients and healthy controls.
- To explore pathway enrichment to understand the functional implications of the previous findings.
- To contextualize findings by comparing results with previous eQTL studies in neurological diseases and reference genotype databases such as the GWAS catalog.

## Repository Structure
The repository is organized as follows:

### `scripts`
Contains all the scripts used in the different analyses.

#### `RNA-seq`
Scripts for RNA-seq analysis to identify DEGs in pDCs from RRMS patients versus healthy controls:
- **`align_run scripts`**: 
  - Performs alignment of RNA-seq data to a reference genome using the `subread-align` tool.
  - Processes paired FastQ files (`_R1` and `_R2`) from a specific directory.
  - Outputs aligned BAM files.
  - Input files (FastQ) are varied for efficiency and to expedite the alignment process.

- **`Read_summarization_featureCounts.R`**: 
  - Utilizes the `Rsubread` library to prepare a list of BAM files.
  - Performs read summarization analysis using `featureCounts`.

- **`RNA_seq_analysis.R`**: 
  - Conducts a comprehensive differential gene expression analysis for plasmacytoid dendritic cells (pDCs) in controls and RRMS patients. Includes:
    - **Quality control**: Library size assessment, removal of low-expression genes, and normalization using TMM.
    - **Exploratory analysis**: Multidimensional scaling (MDS) plots to visualize sample clustering and assess batch effects.
    - **Gene filtering**: Removes genes with low expression or those located on sex chromosomes to reduce noise.
    - **Differential expression analysis**: Uses voom and empirical Bayes statistics to identify DEGs between groups.
    - **Visualization**: Generates heatmaps, MD plots, and other visualizations to represent results.
    - **Exportation**: Exports results (e.g., DEGs) to CSV files for further exploration.

---
