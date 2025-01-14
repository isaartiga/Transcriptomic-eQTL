
---
  title: "Differential gene gene expression analysis of pDCs"
author: "Naiara Garcia Bediaga"
date: "`r Sys.Date()`"
output:
  html_notebook:
  number_sections: yes
toc: yes
toc_float: yes
code_folding: "hide"
html_document:
  toc: yes
number_sections: yes
toc_float: yes
urlcolor: blue
---
  # Introduction
  # This R script performs differential gene expression analysis for plasmacytoid dendritic cells (pDCs)
  # in control individuals (CN) and patients with Relapsing-Remitting Multiple Sclerosis (RRMS).
  # Methods include quality control, normalization, exploratory analysis, and differential expression analysis.
  
  # Load Required Libraries
  # These libraries include core bioinformatics packages for data manipulation, statistical analysis, and visualization.
  

library(data.table)
library(edgeR)
library(limma)
library(BiocGenerics)
library(cowplot)
library(RColorBrewer)
library(gplots)
library(pheatmap)
library(readxl)
library(readr)
library(dplyr)
library(genefilter)
library(grid)
library(biomaRt)
library(fastmatch)
library(ggrepel)
library(fastmatch)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(GO.db)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rtracklayer)
library(KEGG.db)
library(writexl)
library(viridis)
library(compareGroups)


# Load Dataset
# Load the unique read counts per gene, pre-processed and stored in an RDS file.
counts_unique_reads_genes <- readRDS("/mnt/sdb/projects/Iraide/MS_Iraide_Alloza_2023_10/RNA_seq/counts3/Counts_unique_reads_genes_def.rds")

# Load Metadata
# Metadata includes experimental details about each sample, such as group and sequencing run.
samples <- readr::read_csv("/mnt/sdb/projects/Iraide/MS_Iraide_Alloza_2023_10/RNA_seq/RNA-seq/data/metadata_pDC.csv")

# Ensure consistency between metadata and count data
samples <- samples[match(colnames(counts_unique_reads_genes$counts), samples$fileName),]
stopifnot(all(colnames(counts_unique_reads_genes$counts) == samples$fileName))

# Filter only Control and RR samples
filtered_samples <- samples[samples$type %in% c("Control", "RR"), ]
filtered_counts <- counts_unique_reads_genes$counts[, colnames(counts_unique_reads_genes$counts) %in% filtered_samples$fileName]

# Ensure the order matches between counts and metadata
filtered_samples <- filtered_samples[match(colnames(filtered_counts), filtered_samples$fileName), ]

# Quality Control: Library Sizes
# Assess library sizes to ensure uniform sequencing depth across samples.
library_sizes <- data.frame(Library.size = colSums(counts_unique_reads_genes$counts))

# Boxplots of the counts to the power of 0.2 are shown below. The median of all counts is indicated by the dotted line.
# Transform the counts
group <- as.factor(filtered_samples$type)
group.col <- group
levels(group.col) <- brewer.pal(2, "Set2")  # Adjust color palette for two groups (Control and RR)

par(mar = c(8, 4, 4, 2))
z <- filtered_counts^0.2
colnames(z) <- filtered_samples$sampleID  # Use sampleID for the x-axis labels
boxplot(z, main = "Counts^0.2 (Control and RRMS)", col = as.character(group.col), las = 2, ylab = "Transformed Counts")

# Add median line
abline(h = median(filtered_counts^0.2), lty = 2)

# Add legend
legend("topright", 
       legend = levels(group), 
       col = levels(group.col), 
       pch = c(19), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F, 
       inset = c(0.001, 0.001))

# Processing of the data in R
## Annotation of genes
gencode <- as.data.frame(import("/mnt/sdb/references/Genomes/Hg38/Gencode/gencode.v44.primary_assembly.annotation.gtf"))

counts_genes <- counts_unique_reads_genes$counts

gencode1 <- gencode
gencode1 <- gencode1[gencode1$type == "gene", ] 
gencode1 <- gencode1[colSums(is.na(gencode1)) != nrow(gencode1 )]
annotation_genes <- gencode1[match(rownames(counts_genes), gencode1$gene_id),] 
ncbi <- readRDS(file="/mnt/sdb/references/Genomes/Hg38/annotation.library.rds")

annotation_genes<- cbind(annotation_genes, ncbi[match(annotation_genes$gene_name, ncbi$Symbol), c("GeneID", 
                                                                                                  "description", "Synonyms")], Length = annotation_genes$width)

colnames(annotation_genes)[colnames(annotation_genes) == "GeneID"] <- "EntrezID"

# Filter sexual chromosomes and mitochondrial chromosome
# Because the cohort has a mixed of females and males, we are going to filter the genes in the chromosomes X and Y out. We will also filter the chromM out. We only keep Controls and RR samples.
DGE_genes <- DGEList(counts_genes, samples = samples, group = samples$type, genes = annotation_genes)
DGE_genes2 <-  DGE_genes[DGE_genes$genes$seqnames!="chrY",]
DGE_genes2 <-  DGE_genes[DGE_genes$genes$seqnames!="chrY" &  DGE_genes$genes$seqnames!="chrM" &  DGE_genes$genes$seqnames!="chrX",]

# The study is only in RR and controls
DGE_genes2 <- DGE_genes2[,DGE_genes2$samples$type!="PP" ]

## Descriptive analysis of the clinical data
DGE_genes2$samples[DGE_genes2$samples$kit=="miltenyi",]$kit <- "Miltenyi"
DGE_genes2$samples[DGE_genes2$samples$sex=="famale",]$sex <- "female"

# This statistical analysis will help us decide, which variables should be incorporated within the regression models along with the sequencing batch and purification kit. 
export2md(createTable(compareGroups(type~edad+sex+kit+batch,DGE_genes2$samples, min.dis = 1, method = NA, max.ylev = 25, max.xlev=25), extra.labels = c("", "", "", ""),hide.no = c("No",0), show.all=T, show.n = T))


## Exploratory analysis 
### Exploratory analysis of the unadjusted data
# We check the clustering of the samples using a multi-dimensional scaling (MDS) plot. On an MDS plot each sample is plotted according to the distance between it and all other samples. In the MDS plot, the distance between each pair of samples can be interpreted as the leading log-fold change between the samples for the genes that best distinguish that pair of samples. This distance is calculated using the top 500 most variable genes between each pair of samples. The figure above shows that replicate samples from the same group cluster together while samples from different groups are well separated. In other words, differences between groups are much larger than those within groups, meaning that there are likely to be statistically significant differences between the groups. 
# The MDS plots below is showing pre-processed data.They are coloured and shaped according to
# - Colour coded by group and shaped by purification kit
# - Colour coded by group and shaped by sex
# - Colour coded by group and shaped by sequencing run
# - Colour coded by age

plotM=plotMDS(calcNormFactors(DGE_genes2))

temp <- data.frame(Age=DGE_genes2$samples$edad, Group=DGE_genes2$samples$type, kit=DGE_genes2$samples$kit, sex=DGE_genes2$samples$sex,seq_run=DGE_genes2$samples$batch,Leading_logFC_dim_1=plotM$x, Leading_logFC_dim_2=plotM$y)

cat.col <- as.factor(DGE_genes2$samples$type)
levels(cat.col) <- c("#66CDAA", "#EE9572")

ggplot(temp, aes(x=Leading_logFC_dim_1, y=Leading_logFC_dim_2))  +   geom_point(aes(fill=Group, shape=kit), size = 8) +  scale_shape_manual(values = c( 21,22)) + scale_fill_manual(values=levels(cat.col)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) + theme_bw()  + theme(text = element_text(size = 14)) + theme(axis.text = element_text(size = 14))   +  theme(axis.title = element_text(size = 14)) + theme(legend.text = element_text(size = 14)) + xlab(paste0("Leading_logFC_dim_1 (", round(plotM$var.explained[1]*100), "%)")) + ylab(paste0("Leading_logFC_dim_2 (", round(plotM$var.explained[2]*100), "%)")) + ggtitle("MDS coloured by group and shaped by purification kit")

ggplot(temp, aes(x=Leading_logFC_dim_1, y=Leading_logFC_dim_2))  +   geom_point(aes(fill=Group, shape=sex), size = 8) +  scale_shape_manual(values = c( 21,22)) + scale_fill_manual(values=levels(cat.col)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) + theme_bw()  + theme(text = element_text(size = 14)) + theme(axis.text = element_text(size = 14))   +  theme(axis.title = element_text(size = 14)) + theme(legend.text = element_text(size = 14)) + xlab(paste0("Leading_logFC_dim_1 (", round(plotM$var.explained[1]*100), "%)")) + ylab(paste0("Leading_logFC_dim_2 (", round(plotM$var.explained[2]*100), "%)")) + ggtitle("MDS coloured by group and shaped by sex")
# The plot shows that after having removed the sex chromosomes, samples do not cluster by their sex. 

ggplot(temp, aes(x=Leading_logFC_dim_1, y=Leading_logFC_dim_2))  +   geom_point(aes(fill=Group, shape=seq_run), size = 8) +  scale_shape_manual(values = c( 21,22,23,24)) + scale_fill_manual(values=levels(cat.col)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) + theme_bw()  + theme(text = element_text(size = 14)) + theme(axis.text = element_text(size = 14))   +  theme(axis.title = element_text(size = 14)) + theme(legend.text = element_text(size = 14)) + xlab(paste0("Leading_logFC_dim_1 (", round(plotM$var.explained[1]*100), "%)")) + ylab(paste0("Leading_logFC_dim_2 (", round(plotM$var.explained[2]*100), "%)")) + ggtitle("MDS coloured by group and shaped by sequencing run")

ggplot(temp, aes(x=Leading_logFC_dim_1, y=Leading_logFC_dim_2))  +   geom_point(aes(colour = Age), size = 8) +  scale_colour_gradient2(low = "red", mid = "white",high = "blue") + theme_bw()  + theme(text = element_text(size = 14)) + theme(axis.text = element_text(size = 14))   +  theme(axis.title = element_text(size = 14)) + theme(legend.text = element_text(size = 14)) + xlab(paste0("Leading_logFC_dim_1 (", round(plotM$var.explained[1]*100), "%)")) + ylab(paste0("Leading_logFC_dim_2 (", round(plotM$var.explained[2]*100), "%)")) + ggtitle("MDS coloured by age")
# Age does not seem to explain the grouping of the samples. However, age has been associated with gene expression changes and thus we might have to 

### Exploratory analysis of the adjusted data
# MDS plots  with an adjustment included for patient age, and samples´sequencing run and purification kit are shown below. They are coloured according to the group and shapped acoridng to the purification kit. 

logcounts <- cpm(DGE_genes2, log = TRUE)

corrected <- removeBatchEffect(logcounts, batch=DGE_genes2$samples$kit,  batch2=DGE_genes2$samples$batch, covariates =DGE_genes2$samples$edad,  design=model.matrix(~0 + DGE_genes2$samples$group))

colnames(corrected) <- make.unique(DGE_genes2$samples$sampleID, sep = ".")

plotM=plotMDS(corrected)

temp <- data.frame(Age=DGE_genes2$samples$edad ,  Group=DGE_genes2$samples$group, kit=DGE_genes2$samples$kit, sex=DGE_genes2$samples$sex,seq_run=DGE_genes2$samples$batch,Leading_logFC_dim_1=plotM$x, Leading_logFC_dim_2=plotM$y, LB=DGE_genes2$samples$lib.size)

ggplot(temp, aes(x=Leading_logFC_dim_1, y=Leading_logFC_dim_2))  +   geom_point(aes(fill=Group, shape=kit), size = 8) +  scale_shape_manual(values = c( 21,22)) + scale_fill_manual(values=levels(cat.col)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) + theme_bw()  + theme(text = element_text(size = 14)) + theme(axis.text = element_text(size = 14))   +  theme(axis.title = element_text(size = 14)) + theme(legend.text = element_text(size = 14)) + xlab(paste0("Leading_logFC_dim_1 (", round(plotM$var.explained[1]*100), "%)")) + ylab(paste0("Leading_logFC_dim_2 (", round(plotM$var.explained[2]*100), "%)"))

ggplot(temp, aes(x=Leading_logFC_dim_1, y=Leading_logFC_dim_2))  +   geom_point(aes(fill=Group, shape=seq_run), size = 8) +  scale_shape_manual(values = c( 21,22,23,24)) + scale_fill_manual(values=levels(cat.col)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) + theme_bw()  + theme(text = element_text(size = 14)) + theme(axis.text = element_text(size = 14))   +  theme(axis.title = element_text(size = 14)) + theme(legend.text = element_text(size = 14)) + xlab(paste0("Leading_logFC_dim_1 (", round(plotM$var.explained[1]*100), "%)")) + ylab(paste0("Leading_logFC_dim_2 (", round(plotM$var.explained[2]*100), "%)")) + ggtitle("MDS coloured by group and shaped by sequencing run")

ggplot(temp, aes(x=Leading_logFC_dim_1, y=Leading_logFC_dim_2))  +   geom_point(aes(fill=Group), size = 8) +  scale_shape_manual(values = c( 21,22,23,24)) + scale_fill_manual(values=levels(cat.col)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) + theme_bw()  + theme(text = element_text(size = 14)) + theme(axis.text = element_text(size = 14))   +  theme(axis.title = element_text(size = 14)) + theme(legend.text = element_text(size = 14)) + xlab(paste0("Leading_logFC_dim_1 (", round(plotM$var.explained[1]*100), "%)")) + ylab(paste0("Leading_logFC_dim_2 (", round(plotM$var.explained[2]*100), "%)")) + ggtitle("MDS coloured by group and shaped by sequencing run")

## Filtering_genes
# Genes that have very low counts across all the libraries should be removed prior to downstream analysis. This is justified on both biological and statistical grounds. From biological point of view, a gene must be expressed at some minimal level before it is likely to be translated into a protein or to be considered biologically important. From a statistical point of view, genes with consistently low counts are very unlikely be assessed as significantly DE because low counts do not provide enough statistical evidence for a reliable judgement to be made. Such genes can therefore be removed from the analysis without any loss of information.

not.IGgene <- !(grepl("immunoglobulin", DGE_genes2$genes$description) )
DGE_genes2 <- DGE_genes2[not.IGgene, , keep.lib.sizes = FALSE]
nrow(DGE_genes2)

is.exprs_genes <- filterByExpr(DGE_genes2, group = DGE_genes2$samples$group)
DGE_genes2 <- DGE_genes2[is.exprs_genes , , keep.lib.sizes = FALSE]
table(is.exprs_genes)

## Normalization_genes with TMM 
# Normalization by trimmed mean of M values (TMM) (Robinson and Oshlack 2010) is performed by using the calcNormFactors function, which returns the DGEList argument with only the norm.factors changed.The scaling factor for each sample:

DGE_genes2 <- normLibSizes(DGE_genes2)
logcpm <- cpm(DGE_genes2, log = TRUE)

write.csv(logcpm, "logcounts_pDCs_RRvsHC.csv")

## Biological variation 
# The biological variation between replicate samples is estimated, the common BCV:
design <- model.matrix(~0 + DGE_genes2$samples$group + DGE_genes2$samples$kit + DGE_genes2$samples$batch + DGE_genes2$samples$edad)
colnames(design ) <- c(levels(as.factor(DGE_genes2$samples$group)),"StemCell", levels(as.factor(DGE_genes2$samples$batch))[-1], "edad")

disp_genes <- estimateDisp(DGE_genes2, design = design, robust = TRUE)
sqrt(disp_genes$common.dispersion)

# A summary of the gene specific dispersions:
summary(sqrt(disp_genes$tagwise.dispersion))


# The BCV plot with highly variable genes highlighted:
#  The vertical axis of the plotBCV plot shows square-root dispersion, also known as biological coefficient of variation (BCV) (McCarthy, Chen, and Smyth 2012). For RNA-seq studies, the NB dispersions tend to be higher for genes with very low counts. 
#- *Common dispersion* is the mean dispersion across all genes.
#- *Trended dispersion* is the mean dispersion across all genes with similar abundance. In other words, the fitted value of the mean-dispersion trend.
# Now, as for how this affects hypothesis testing - there are three main points, in order of obviousness:
# -Larger dispersions = higher variance between replicates, which reduce power to detect DE.
# - The performance of the model (and thus of the DE analysis) depends on the accuracy of the dispersion estimate. If there is a strong mean-dispersion trend, the common dispersion is obviously unsuitable. If the gene-specific dispersions vary around the trend, the trended dispersion is unsuitable. The "raw" tagwise estimates are unbiased estimates of the gene-specific dispersions, and should be the most suitable, except...
# The performance of the model also depends on the precision with which the dispersions are estimated. Here, the raw estimates are least stable as they use the least amount of information, whereas the trended (and to a greater extent, common) dispersion estimates share information between genes for greater stability. This is why the shrunken tagwise estimates (that you get with default estimateDisp) are so useful, as they provide a compromise between precision and accuracy.

plotBCV(disp_genes)

HV_genes <- disp_genes[disp_genes$prior.df == min(disp_genes$prior.df), ]


## Differential Gene Expression Analysis
# This section evaluates differential expression between the experimental groups (RR vs Control).

# Defining the contrast:
# Specify the contrast of interest between the Relapsing-Remitting (RR) and Control groups.
cont_m <- makeContrasts(RRvsCtr = RR - Control, levels = design)

# Data transformation with voom:
# Gene expression data is transformed to log2-CPM, and precision weights are calculated for each observation.
# These weights adjust the expected variance according to the mean, improving model robustness.
vQW <- voomWithQualityWeights(DGE_genes2, design = design, plot = TRUE)

# Linear model fitting:
# A linear model is fitted for each gene using the voom-transformed data.
fit_genes <- lmFit(vQW, design = design)

# Applying the contrast:
# The previously defined contrast is applied to the fitted model to compare groups.
fit_genes <- contrasts.fit(fit_genes, contrasts = cont_m)

# Bayesian statistics:
# Empirical Bayes moderation is applied to stabilize the variance estimates.
fit_genes <- eBayes(fit_genes)

# Plotting the mean-variance trend:
# This plot visualizes the mean-variance relationship after model adjustments.
plotSA(fit_genes, main = "Final model: Mean-variance trend")

# Deciding significant tests:
# Identify significantly differentially expressed genes with a p-value threshold of 0.05.
DT <- decideTests(fit_genes, p.value = 0.05)
summary(DT)

# RR vs Control visualization:
# MD plot: Visualizes log2 fold-change vs average expression for the RR vs Control comparison.
plotMD(fit_genes, coef = 1, status = DT[, "RRvsCtr"], hl.col = c("red", "blue"), 
       values = c(1, -1), legend = FALSE, xlab = "Average log CPM", 
       ylab = "Log-fold-change", cex.lab = 1.5)
legend("topright", c("Up", "Down", "NotSig"), pch = 16, 
       col = c("red", "blue", "black"), pt.cex = c(1, 1, 0.2))

# Heatmap visualization:
# Create a heatmap for the top differentially expressed genes.

# Calculate log2-transformed counts.
logcounts_filt <- cpm(DGE_genes2, log = TRUE)

# Remove batch effects to focus on biological variation.
corrected_filt <- removeBatchEffect(logcounts_filt, batch2 = DGE_genes2$samples$batch, 
                                    covariates = DGE_genes2$samples$edad, 
                                    design = model.matrix(~0 + DGE_genes2$samples$group))

# Annotate rows and columns for the heatmap.
row.names(corrected_filt) <- DGE_genes2$genes$gene_name
colnames(corrected_filt) <- make.unique(DGE_genes2$samples$sampleID)

# Define color palette and group annotations for the heatmap.
HM_colorp <- colorRampPalette(c("blue3", "white", "red3"))
HM_col_anno <- data.frame(Group = as.factor(DGE_genes2$samples$group), 
                          row.names = colnames(corrected_filt))

group.col <- as.factor(DGE_genes2$samples$group)
levels(group.col) <- brewer.pal(3, "Set2")[-3]
HM_cols <- list(Group = levels(group.col))
names(HM_cols$Group) <- levels(as.factor(DGE_genes2$samples$group))

# Generate the heatmap for the top 100 differentially expressed genes.
pheatmap(corrected_filt[row.names(corrected_filt) %in% tt_RRvsCtr$gene_name[1:100], ],
         color = HM_colorp(256), 
         scale = "row", 
         cluster_cols = TRUE,  
         annotation_col = HM_col_anno, 
         annotation_colors = HM_cols,
         fontsize_row = 6)

# Extract the top differentially expressed genes:
# Filter genes with adjusted p-value < 0.05 and save the results to a CSV file.
tt_RRvsCtr <- topTable(fit_genes, coef = "RRvsCtr", n = Inf, sort.by = "p")
DEG_bulk <- tt_RRvsCtr[tt_RRvsCtr$adj.P.Val < 0.05,]
write_csv(as.data.frame(DEG_bulk$gene_name), "DEG_bulk_pDCs_RRvsHC.csv")

