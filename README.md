# Transcriptomic-eQTL

## Overview
This repository contains the scripts and results for an eQTL analysis conducted as part of a Master's Thesis. The study investigates the regulatory effects of genetic variants on gene expression in a cohort of Relapsing-Remitting Multiple Sclerosis (RRMS) patients and healthy controls (HC).

## Objectives
- **Differential Gene Expression Analysis**: Identify transcriptional differences between RRMS patients and healthy controls (HC) by discovering differentially expressed genes (DEGs).
- **Genotyping Quality Control (QC)**: Perform rigorous QC procedures on genotyping data to ensure the accuracy and reliability of genetic variants and sample integrity for cis-eQTL analyses.
- **cis-eQTL Analysis**: Uncover genetic variants influencing gene expression in plasmacytoid dendritic cells (pDCs) using DEGs.
- **Genotypic Variability**: Compare the allelic frequencies of key cis-eQTL-associated single nucleotide polymorphisms (SNPs) between RRMS patients and healthy controls.
- **Pathway Enrichment**: Explore the functional implications of cis-eQTL results through GO and KEGG pathway enrichment.

## Repository Structure
The repository is organized into folders and scripts according to the analysis workflows:

### **Workflow Diagram**
A comprehensive workflow diagram (found in `workflow_transcriptomic_eQTL.pdf`) summarizes the analyses, bioinformatic tools used, and results.

### **Folders**
- **`scripts/`**: Contains R and shell scripts for RNA-Seq, genotype QC, eQTL analysis, and post-eQTL analysis.
- **`results/`**: Includes the DEG and significant cis-eQTL tables, along with diagnostic plots and pathway enrichment results.

## Key Workflows

### RNA-Seq Analysis
This analysis identifies DEGs in pDCs between RRMS and HC samples:
- **Alignment**: Scripts align RNA-seq data to a reference genome using `subread-align`.
- **Read Summarization**: Generates gene count tables using the `featureCounts` function in `Rsubread`.
- **Differential Expression**: Conducts DEG analysis with the following steps:
  - **Quality Control**: Normalization and removal of low-expression genes.
  - **Exploratory Analysis**: MDS plots for clustering and batch effect evaluation.
  - **DEG Analysis**: Using `voom` and empirical Bayes methods.
  - **Visualization**: Heatmaps, MD plots, and other diagnostic plots.

*Note: The RNA-Seq analysis scripts were initially created by the bioinformatics platform at IIS Biobizkaia. These scripts have been adapted for this study and tailored to meet the specific needs of this research.*

### Genotype Quality Control (QC)
The **`genotype_QC.sh`** pipeline performs quality control (QC) at both the variant and sample levels. It calls individual R scripts throughout the process to:
- Filter variants by:
  - **Missingness**: Excludes SNPs with missingness rates > 0.05.
  - **Case-control rate differences**: Removes SNPs with significant differences (p < 1 × 10⁻⁵).
  - **Minor Allele Frequency (MAF)**: Filters SNPs with MAF < 0.05 to reduce false-positive findings.
  - **Hardy-Weinberg Equilibrium (HWE)**: Excludes SNPs with HWE violations (p < 1 × 10⁻¹⁰ in cases, p < 1 × 10⁻⁶ in controls).
  - **Chromosomal Location**: Excludes SNPs outside chromosomes 1–23.

- Filter samples by:
  - **Missingness**: Removes samples with missingness rates > 0.05.
  - **Heterozygosity**: Detects and removes samples with deviations in heterozygosity rates.
  - **Sex Discrepancies**: Identifies mismatches between reported and genetic sex.
  - **Relatedness**: Excludes samples with PI_HAT > 0.2.
  - **Population Stratification**: Uses MDS analysis to identify and exclude population outliers by merging with the 1000 Genomes Project reference dataset.

*Note: The genotype QC scripts were also developed by the bioinformatics platform at IIS Biobizkaia and customized for this analysis.*

### cis-eQTL Analysis
A 1Mb window is applied around DEGs for cis-eQTL discovery using the `MatrixeQTL` package:
- **Data Preparation**: Scripts normalize expression data, format genotype files, and filter SNPs within 1Mb of DEGs.
- **Analysis**: Uses the `modelLINEAR` model to identify significant cis-eQTLs (FDR < 0.05).
- **Annotation and Visualization**: Annotates results and visualizes genotype-expression associations with violin plots.

### Post-eQTL Analysis
Explores the biological significance of cis-eQTL findings:
- **Genotypic Variability**: Compares allelic frequencies between RRMS and HC.
- **Functional Enrichment**: Conducts GO and KEGG pathway analysis to interpret functional relevance.

## Results
The `results/` folder contains:
- **DEG Results**: Complementary Table 1 with complete DEG information.
- **cis-eQTL Results**: Complementary Table 2 with significant cis-eQTL findings.

## Usage
1. **RNA-Seq Analysis**:
   - Run the scripts sequentially: `align_run`, `Read_summarization_featureCounts.R`, and `RNA_seq_analysis.R`.
2. **Genotype QC**:
   - Execute `genotype_QC.sh` to process and filter genotyping data.
3. **cis-eQTL Analysis**:
   - Prepare data using the formatting scripts, then run `run_matrix_eQTL_final.R`.
4. **Post-eQTL Analysis**:
   - Use scripts in `post_eQTL_analysis` for variability and enrichment studies.

