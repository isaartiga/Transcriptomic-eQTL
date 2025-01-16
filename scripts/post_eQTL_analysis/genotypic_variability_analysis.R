# Load genotype data and cis-eQTL results) 
genotype <- read.table("SNP_DEG.txt", sep = "\t", header = TRUE)
cis_FDR_pDCs <- read.csv("cis_FDR_pDCs.csv")

# Filter SNPs based on those present in the cis_FDR_pDCs dataset
genotype_pDCs <- genotype[genotype$vcf.ID %in% cis_FDR_pDCs$snps, ]

# Load necessary libraries
library(dplyr)

# Identifying patient (RR) and control (CN) columns in the dataset
patients <- grep("^RR", colnames(genotype_pDCs), value = TRUE)
controls <- grep("^CN", colnames(genotype_pDCs), value = TRUE)

# Create a dataframe to store statistical test results
test_RRvsCN_genotype_pDCs <- data.frame(SNP = character(),
                                     Test = character(),
                                     P_value = numeric(),
                                     stringsAsFactors = FALSE)

# Loop through each SNP (row) in the dataset
for (snp in 1:nrow(genotype_pDCs)) {
    # Extract SNP name
    snp_name <- genotype_pDCs[snp, "vcf.ID"]
    
    # Count occurrences of 0, 1, 2 for patients and controls
    count_patients <- table(factor(as.numeric(genotype_pDCs[snp, patients]), levels = c(0, 1, 2)))
    count_controls <- table(factor(as.numeric(genotype_pDCs[snp, controls]), levels = c(0, 1, 2)))
    
    # Create a contingency table for the genotype counts
    contingency_table <- rbind(count_patients, count_controls)
    
    # Decide the statistical test to use
    if (any(contingency_table < 5)) {
        test <- fisher.test(contingency_table)  # Use Fisher's exact test if expected counts are low
        test_type <- "Fisher"
    } else {
        test <- chisq.test(contingency_table)  # Use Chi-squared test otherwise
        test_type <- "Chi-squared"
    }
    
    # Store the test results in the dataframe
    test_RRvsCN_genotype_pDCs <- rbind(test_RRvsCN_genotype_pDCs, 
                                    data.frame(SNP = snp_name,
                                               Test = test_type,
                                               P_value = test$p.value))
}

# Filter results for significant p-values (p < 0.05)
test_RRvsCN_genotype_pDCs_pval <- test_RRvsCN_genotype_pDCs[test_RRvsCN_genotype_pDCs$P_value < 0.05,]

# Handle SNP-to-multiple-gene associations and sort them alphabetically
test_RRvsCN_genotype_pDCs_pval$Gene <- sapply(
    test_RRvsCN_genotype_pDCs_pval$SNP, 
    function(snp) {
        # Filter genes associated with the SNP in the cis_FDR_pDCs_annotated dataset
        associated_genes <- cis_FDR_pDCs_annotated$Gene[cis_FDR_pDCs_annotated$rs == snp]
        
        # Sort genes alphabetically and concatenate them with commas
        if (length(associated_genes) > 0) {
            return(paste(sort(associated_genes), collapse = ", "))
        } else {
            return(NA)  # Return NA if no associated genes are found
        }
    }
)

# Function to calculate the risk allele frequency for a given SNP
Risk_Allele_Frequency <- function(row, samples) {
    alleles <- as.numeric(genotype_pDCs[row, samples])
    freq <- sum(alleles, na.rm = TRUE) / (2 * length(samples)) * 100  # Calculate allele frequency
    return(freq)
}

# Calculate risk allele frequencies for patients and controls
Risk_Allele_Frequency_RRMS <- sapply(1:nrow(genotype_pDCs), Risk_Allele_Frequency, samples = patients)
Risk_Allele_Frequency_HC <- sapply(1:nrow(genotype_pDCs), Risk_Allele_Frequency, samples = controls)

# Combine allele frequencies and statistical results into a final dataframe
final_df_dif_genotype <- data.frame(
    `SNP` = test_RRvsCN_genotype_pDCs_pval$SNP,
    `Gene` = test_RRvsCN_genotype_pDCs_pval$Gene,
    `p-value` = test_RRvsCN_genotype_pDCs_pval$P_value,
    `Risk Allele Frequency (RRMS)` = Risk_Allele_Frequency_RRMS[match(test_RRvsCN_genotype_pDCs_pval$SNP, genotype_pDCs$vcf.ID)],
    `Risk Allele Frequency (HC)` = Risk_Allele_Frequency_HC[match(test_RRvsCN_genotype_pDCs_pval$SNP, genotype_pDCs$vcf.ID)]
)

# Add the risk allele for each SNP
final_df_dif_genotype$Risk_Allele <- cis_FDR_pDCs_annotated$Alt.allele[match(final_df_dif_genotype$SNP, cis_FDR_pDCs_annotated$rs)]

# Save the results to an Excel file
library(writexl)
write_xlsx(final_df_dif_genotype, "results_tables/tabla_final_dif_genotype.xlsx")
