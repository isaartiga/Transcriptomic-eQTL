# Load necessary libraries for reading CSV and Excel files
library(readr)
library(readxl)

# Load logcounts data from a CSV file and samplesheet data for metadata
logcounts <- read_csv("pDC_logCounts_Isabel.csv")
samplesheet <- read_csv("conversion_samplesheet_pDCs_64.csv")

# Step 1: Create a mapping between fileName and sampleID in the samplesheet
# This mapping will allow renaming logcounts columns from fileName to sampleID
mapping <- setNames(samplesheet$sampleID, samplesheet$fileName)

# Rename the columns in logcounts (excluding the first column of gene names)
# using the mapping created above
colnames(logcounts)[-1] <- mapping[colnames(logcounts)[-1]]

# Step 2: Ensure `samplesheet` has the columns `fileName` and `sampleID`
# Create a vector of sample names in the desired order based on `sampleID` in `samplesheet`
desired_order <- samplesheet$sampleID

# Step 3: Reorder the columns of logcounts to match the order of `sampleID` in `samplesheet`
logcounts <- logcounts[, c(1, match(desired_order, colnames(logcounts)[-1]) + 1)]

# Extract necessary columns from the samplesheet for batch correction
kit <- samplesheet$kit                # Variable representing the first batch effect
set <- samplesheet$batch              # Variable representing the second batch effect (if any)
edad <- samplesheet$edad              # Covariate for age
group <- factor(samplesheet$type)     # Experimental group, converted to a factor for design matrix

# Ensure that the first column of logcounts contains gene names
genes <- logcounts[, 1]               # Extract gene names
expr_data <- logcounts[, -1]          # Remove the first column to get expression data only

# Create the design matrix based on the group factor for batch correction
design <- model.matrix(~0 + group)

# Step 4: Install the `limma` package if it is not already installed
# This package contains the function `removeBatchEffect` for batch correction
if (!requireNamespace("limma", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install("limma")
}

# Load the `limma` package
library(limma)

# Step 5: Apply batch correction using removeBatchEffect
# We correct for two batch variables (kit and set) and the experimental group design
expr_corrected <- removeBatchEffect(
    expr_data,
    batch = kit,
    batch2 = set,
    design = design
)

# Load DEG (differentially expressed genes) data from an Excel file
library(readxl)
DEG <- read_excel("todos_DEG_BULK.xlsx")

# Step 6: Combine corrected expression data with gene names in a data frame
logcounts_corrected <- as.data.frame(cbind(genes, expr_corrected))
colnames(logcounts_corrected)[1] <- "Gene"  # Rename the first column to "Gene"

# Filter logcounts_corrected to only include genes present in the DEG list
logcounts_corrected_DEG <- logcounts_corrected[logcounts_corrected$Gene %in% DEG$gene_name,]

# Save the filtered corrected expression data to a TSV file
write_tsv(logcounts_corrected_DEG, "GE_corrected.txt")

# Step 7: Load gene annotation data for chromosome and location information
genes_annot <- read_csv("pDC_genes_annot.csv")

# Step 8: Extract chromosome and location information based on gene names in DEG data
# Map gene IDs in logcounts_corrected_DEG to the annotation data in genes_annot
geneid <- logcounts_corrected_DEG$Gene
chr <- genes_annot$seqnames[match(geneid, genes_annot$gene_name)]
right <- genes_annot$start[match(geneid, genes_annot$gene_name)]
left <- genes_annot$end[match(geneid, genes_annot$gene_name)]

# Step 9: Create a data frame containing gene location information
geneloc <- as.data.frame(cbind(geneid, chr, right, left))

# Save gene location data to a TSV file
write_tsv(geneloc, "geneloc.txt")




