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

- *scripts* : todos los scripts de los distintos análsisi estan en esta carpeta
  - *RNA-seq* : contiene los scripts del analisis de RNA-seq para determinar DEGs en pDCs de RRMS vs HC. 
    - *align_run* scripts:  realiza la alineación de datos de secuenciación de ARN (RNA-seq) a un genoma de referencia utilizando la herramienta subread-align, procesando archivos FastQ emparejados (_R1 y _R2) en un directorio específico, y generando archivos BAM alineados como salida.únicamente cambian los archivos de entrada (fastq) para agilizar el proceso.
    - *Read_summarization_featureCounts.R* : utiliza la librería Rsubread para preparar una lista de archivos BAM para un análisis posterior de sumarización de lecturas (read summarization) con featureCounts.
    - *RNA_seq_analysis.R* : realiza un análisis de expresión diferencial para células dendríticas plasmacitoides (pDCs) en individuos de control (Control) y pacientes con esclerosis múltiple remitente-recurrente (RRMS). Incluye control de calidad, análisis exploratorio de variables, filtrado y normalización de genes, analisis de expresión diferencial, incluyendo exportación y visualización de resultados. 
