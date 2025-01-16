# Transcriptomic-eQTL

## Overview
This repository contains the scripts and results for an eQTL analysis conducted as part of a Master's Thesis. The study investigates the regulatory effects of genetic variants on gene expression in a cohort of Relapsing-Remitting Multiple Sclerosis (RRMS) patients and healthy controls (HC).

## Objectives
- **Differential Gene Expression Analysis**: Identify transcriptional differences between RRMS patients and healthy controls (HC) by discovering differentially expressed genes (DEGs).
- **Genotyping Quality Control (QC)**: Perform rigorous QC procedures on genotyping data to ensure the accuracy and reliability of genetic variants and sample integrity for cis-eQTL analyses.
- **cis-eQTL Analysis**: Uncover genetic variants influencing gene expression in plasmacytoid dendritic cells (pDCs) using DEGs.
- **Genotypic Variability**: Compare the allelic frequencies of key cis-eQTL-associated single nucleotide polymorphisms (SNPs) between RRMS patients and healthy controls.
- **Pathway Enrichment**: Explore the functional implications of cis-eQTL results through GO and KEGG pathway enrichment.

A comprehensive workflow diagram summarizes the analyses, bioinformatic tools used, and results (workflow_transcriptomic_eQTL).


---

## Repository Structure
All scripts are located in the `scripts` folder and are categorized by analysis type: RNA-Seq, genotype_QC, eQTL_analysis, and post_eQTL_analysis (Genotypic Variability and Functional Enrichment). Below is an overview of the structure and key workflows:

### RNA-Seq Analysis
This analysis identifies DEGs in pDCs between RRMS and HC samples using the following scripts:
- **`align_run scripts`**:
  - Aligns RNA-seq data to a reference genome using the `subread-align` tool.
  - Processes paired FastQ files (`_R1` and `_R2`) for each sequencing run and generates BAM files.
- **`Read_summarization_featureCounts.R`**:
  - Summarizes read counts from BAM files using the `Rsubread` library and prepares data for DEG analysis.
- **`RNA_seq_analysis.R`**:
  - Performs differential expression analysis with steps including:
    - **Quality Control**: Library size checks, removal of low-expression genes, normalization (TMM).
    - **Exploratory Analysis**: MDS plots for clustering and batch effect evaluation.
    - **Gene Filtering**: Removes genes on sex chromosomes.
    - **DEG Analysis**: Utilizes `voom` and empirical Bayes methods.
    - **Visualization**: Generates heatmaps, MD plots, and other visualizations.
    - **Export**: Outputs DEGs as CSV files.

---

### Genotype Quality Control (QC)
A main pipeline script (**`genotype_QC.sh`**) performs quality control (QC) at both the variant and sample levels. It calls individual R scripts throughout the process.

---

### cis-eQTL Analysis
A 1Mb window is applied around DEGs for cis-eQTL discovery using the `MatrixeQTL` package. Input files are formatted specifically for the analysis:
- **`prepare_expression_data.R`**: Normalizes expression data and corrects for batch effects (sequencing kit and run).
- **`vcf2matrixeqtlformat.R`**: Converts VCF files from genotype QC into the required format.
- **`filter_SNP_GE.R`**: Filters SNPs and genes within 1Mb windows around DEGs.
- **`prepare_covariates_files.R`**: Creates a covariates file using age as a covariate.
- **`run_matrix_eQTL_final.R`**:
  - Conducts cis-eQTL analysis using the `modelLINEAR` model.
  - Filters significant cis-eQTLs by FDR < 0.05.
- **`annotation_cis_eQTL.R`**: Annotates significant cis-eQTL results with genomic features.
- **`graficos_violin_eqtl.R`**: Creates violin plots visualizing expression levels versus genotypes.

---

### Post-eQTL Analysis
This analysis explores the biological significance of cis-eQTL findings:
- **`genotypic_variability_analysis.R`**:
  - Compares allelic frequencies between RRMS and HC using Fisher’s test.
- **`GRAFICO_GENOTIPO_EXP_RR_CN.R`**:
  - Visualizes genotype-expression associations using scatter plots.
- **`script_go_kegg_final.R`**:
  - Performs functional enrichment analysis (GO and KEGG).
  - Outputs pathway enrichment plots.

---

### Results
The `results` folder contains the complete DEG and significant cis-eQTL tables, excluded from the main thesis due to length.

---

## Usage
1. **RNA-Seq Analysis**:
   - Run the `align_run`, `Read_summarization_featureCounts.R`, and `RNA_seq_analysis.R` scripts sequentially.
2. **Genotype QC**:
   - Execute `genotype_QC.sh` to process and filter genotyping data.
3. **cis-eQTL Analysis**:
   - Prepare data using the formatting scripts, then run `run_matrix_eQTL_final.R`.
4. **Post-eQTL Analysis**:
   - Explore variability and enrichment using scripts in the `post_eQTL_analysis` folder.

---



