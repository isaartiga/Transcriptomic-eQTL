# Load gene location data from the geneloc file
geneloc <- read.table("geneloc.txt", header = TRUE, sep = "\t")

# Load dplyr library for data manipulation
library(dplyr)

# Convert 'right' and 'left' columns to numeric (if they're read as non-numeric data types)
geneloc$right <- as.numeric(as.character(geneloc$right))
geneloc$left <- as.numeric(as.character(geneloc$left))

# Step 1: Define a range of ±1 million base pairs around each gene
# Adds two new columns, `start_range` and `end_range`, that define a window around each gene
geneloc_window <- geneloc %>%
    mutate(start_range = right - 1e6,
           end_range = left + 1e6)

# Load SNP data from SNP_all.txt and rename the first column to "vcf.ID" for consistency
SNP <- read.table("SNP_all.txt", header = TRUE, sep = "\t")
colnames(SNP)[1] <- "vcf.ID"

# Load SNP position data from SNPpos_all.txt
snpspos <- read.table("SNPpos_all.txt", header = TRUE, sep = "\t")

# Step 2: Filter SNPs within ±1 million base pairs of genes in `geneloc_window`
# Join `snpspos` with `geneloc_window` on chromosome (`chr`) and filter by position range
# to keep SNPs within the defined range. Remove duplicates based on `snpid`.
filtered_variants_DEG <- snpspos %>%
    inner_join(geneloc_window, by = "chr", relationship = "many-to-many") %>%
    filter(pos >= start_range & pos <= end_range) %>%
    distinct(snpid, .keep_all = TRUE) %>%  # Remove duplicates by `snpid`
    select(snpid, chr, pos)  # Keep only columns of interest

# Step 3: Filter the main SNP dataset to include only SNPs found in `filtered_variants_DEG`
# Select SNPs in `SNP` with `vcf.ID` matching `snpid` in `filtered_variants_DEG`
SNP_bulk <- SNP[SNP$vcf.ID %in% filtered_variants_DEG$snpid, ]

# Remove duplicate SNPs in `SNP_bulk`, ensuring unique `vcf.ID` rows
SNP_bulk <- SNP_bulk %>%
    distinct(vcf.ID, .keep_all = TRUE)

# Save the filtered SNP data to a TSV file
library(readr)
write_tsv(SNP_bulk, "SNP_DEG.txt")

# Step 4: Create and save a subset of SNP position data for SNPs in SNP_bulk
# Filter `snpspos` to include only SNPs present in `SNP_bulk` and save it
snpspos_bulk <- snpspos[snpspos$snpid %in% SNP_bulk$vcf.ID, ]
write_tsv(snpspos_bulk, "snpsloc_DEG.txt")

