# Load metadata and expression
GE <- read.table("GE_corrected.txt", sep = "\t", header = TRUE)
samplesheet <- read.csv("conversion_samplesheet_pDCs_64.csv")

# Extract sample IDs from the column names of the gene expression data (excluding "Gene" column)
sampleid <- colnames(GE)[-1]

# Match each sample ID to its age using the samplesheet, based on `sampleID`
age <- samplesheet$edad[match(sampleid, samplesheet$sampleID)]

# Convert age data to a data frame, transposing so it matches the format needed for covariates
covariates_age <- as.data.frame(t(age))
colnames(covariates_age) <- sampleid  # Set column names as sample IDs

# Add an ID row label for "age" to indicate this row represents age data
covariates_age <- cbind("ID" = "age", covariates_age)

# Save the age covariates data to a file
library(readr)
write_tsv(covariates_age, "cov_age.txt")
