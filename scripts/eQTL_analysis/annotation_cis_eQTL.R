# This script annotates SNPs with relevant gene information for significative cis-eQTL. 
# It merges SNP data with gene annotations, calculates the distance between each SNP and its associated gene, 
# and prepares a cleaned and formatted dataset for downstream analysis.

# Load the gene annotation data
genes_annot <- read_csv("/mnt/sdb/iartiga/analisis_MatrixeQTL_final/pDC_genes_annot.csv")

# Merge SNP data (`axiom`) with eQTL results (`cis_FDR_ESTRATEGIA_A`) based on SNP IDs
# Also merge with gene annotation data to include relevant gene information
cis_FDR_A_anotado <- merge(axiom, cis_FDR_ESTRATEGIA_A, by.y = "snps", by.x = "Extended RSID", all.x = FALSE)
cis_FDR_A_anotado <- merge(cis_FDR_A_anotado, genes_annot[, c(2, 3, 8, 10)], by.x = "gene", by.y = "gene_name")

# Define a function to calculate the distance between a SNP and a gene
# The function takes the SNP position, gene start position, and gene end position as inputs
calculate_distance <- function(snp_pos, gen_start, gen_end) {
    if (snp_pos < gen_start) {
        return(gen_start - snp_pos)  # SNP is located before the gene
    } else if (snp_pos > gen_end) {
        return(snp_pos - gen_end)    # SNP is located after the gene
    } else {
        return(0)  # SNP is located within the gene
    }
}

# Ensure positional columns are numeric for calculations
cis_FDR_A_anotado$`Physical Position` <- as.numeric(as.character(cis_FDR_A_anotado$`Physical Position`))
cis_FDR_A_anotado$start <- as.numeric(as.character(cis_FDR_A_anotado$start))
cis_FDR_A_anotado$end <- as.numeric(as.character(cis_FDR_A_anotado$end))

# Calculate the distance between each SNP and its associated gene
# The `mapply` function applies `calculate_distance` to each row of the dataframe
cis_FDR_A_anotado$Distance <- mapply(
    calculate_distance, 
    cis_FDR_A_anotado$`Physical Position`, 
    cis_FDR_A_anotado$start, 
    cis_FDR_A_anotado$end
)

# Select and reorder relevant columns for the final annotated dataframe
cis_FDR_A_anotado <- cis_FDR_A_anotado[, c(
    "Extended RSID", "Chromosome", "Physical Position", "Ref Allele", "Alt Allele", 
    "Distance", "gene", "gene_id", "start", "end", "statistic", "pvalue", "FDR", "beta"
)]

# Rename column names to more user-friendly labels
colnames(cis_FDR_A_anotado) <- c(
    "rs", "Chromosome", "SNP.pos", "Ref.allele", "Alt.allele", 
    "Distance2gene", "Gene", "Gene_ID", "Gene.start", "Gene.end", 
    "Statistic", "pvalue", "FDR", "beta"
)

# Order the annotated data by FDR values in ascending order
cis_FDR_A_anotado <- cis_FDR_A_anotado[order(cis_FDR_A_anotado$FDR),]

# Load the dplyr library for additional data manipulation
library(dplyr)

# Remove duplicate entries based on the `rs` (SNP ID) and `Gene` columns
cis_FDR_A_anotado <- cis_FDR_A_anotado %>%
    distinct(rs, Gene, .keep_all = TRUE)
