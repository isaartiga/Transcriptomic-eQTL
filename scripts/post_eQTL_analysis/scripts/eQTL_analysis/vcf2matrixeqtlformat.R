# Load SNP data from a VCF file, ignoring comment lines starting with "#"
SNP <- read.table("ALL_passed_QC.vcf", comment.char = "#", header = FALSE, sep = "\t")

# Read all lines of the VCF file to capture column names
SNP_lines <- readLines("ALL_passed_QC.vcf")

# Find the line that defines the column names, starting with "#CHROM"
line_colnames <- grep("^#CHROM", SNP_lines, value = TRUE)

# Remove the leading "#" from the column names line
line_colnames <- sub("^#", "", line_colnames)

# Split the line to extract the column names and assign them to SNP
colnames(SNP) <- strsplit(line_colnames, "\t")[[1]]

# Extract SNP IDs and keep only columns with genotypes (starting from column 10)
snpid <- SNP$ID
SNP <- cbind(snpid, SNP[, 10:ncol(SNP)])

# Recode genotypes to numeric format
SNP <- as.data.frame(lapply(SNP, function(x) {
  x <- gsub("1/1", "2", x)   # Convert "1/1" to "2"
  x <- gsub("0/1", "1", x)   # Convert "0/1" to "1"
  x <- gsub("0/0", "0", x)   # Convert "0/0" to "0"
  x <- gsub("1/0", "1", x)   # Convert "1/0" to "1"
  x <- gsub("./.", NA, x)    # Convert missing values "./." to NA
  return(x)
}))

# Convert SNP columns from character to numeric, handling "." as NA
SNP[-1] <- lapply(SNP[, 2:ncol(SNP)], function(x) {
  if (is.character(x)) {
    x[x == "."] <- NA
    as.numeric(x)
  } else {
    x
  }
})

# Clean column names by transforming them to a consistent format
colnames_SNP_clean <- gsub("([A-Za-z]+)(\\d+)[_\\.](\\d+).*", "\\1\\2-\\3", colnames(SNP))
colnames(SNP) <- colnames_SNP_clean

# Load the samplesheet containing the mapping of file names to sample IDs
samplesheet <- read_csv("/mnt/sdb/iartiga/analisis_MatrixeQTL_final/conversion_samplesheet_pDCs_64.csv")

# Create a mapping from KV (file names) to sampleID
mapping <- setNames(samplesheet$sampleID, samplesheet$KV)

# Rename the columns in SNP (excluding the first column with SNP IDs) based on this mapping
colnames(SNP)[-1] <- mapping[colnames(SNP)[-1]]

# Define the desired order of samples based on sampleID in the samplesheet
desired_order <- samplesheet$sampleID

# Filter `desired_order` to include only samples present in `SNP`
filtered_order <- desired_order[desired_order %in% colnames(SNP)[-1]]

# Reorder SNP to match the filtered sample order
SNP <- SNP[, c(1, match(filtered_order, colnames(SNP)[-1]) + 1)]

# Save the processed SNP data to a TSV file
library(readr)
write_tsv(SNP, "SNP_all.txt")

# Create SNP location data by extracting chromosome and position information from VCF
SNP_loc <- as.data.frame(cbind(snpid, SNP$CHROM, SNP$POS))

# Rename columns in SNP_loc to "chr" and "pos" for chromosome and position information
colnames(SNP_loc) <- c("snpid", "chr", "pos")

# Prefix chromosome numbers with "chr" for consistency
SNP_loc$chr <- paste0("chr", SNP_loc$chr)

# Save the SNP location data to a TSV file
write_tsv(SNP_loc, "SNPpos_all.txt")
