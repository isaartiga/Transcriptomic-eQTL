# Common configuration
# This script performs an eQTL analysis using the MatrixEQTL package, specifically focusing on the 'STRATEGY_A'.
# It is designed to be modular and scalable for future strategies.

# Install MatrixEQTL if not already installed
if (!require(MatrixEQTL)) {
  install.packages("MatrixEQTL", repos = "http://cran.r-project.org")
  library(MatrixEQTL)
} else {
  library(MatrixEQTL)  # Load MatrixEQTL library for eQTL analysis
}

# Base directory where input and output files are stored
base_directory <- "/mnt/sdb/iartiga/FINAL_ANALYSIS_eQTL/"

# Configuration for pDCs
# This list stores the parameters for the analysis, including the gene expression file and covariates.
analysis_config <- list(
  list(name = "pDCs", gene_expression_file = "GE_corrected.txt", covariates_file = "cov_age.txt")
)

# Initialize an empty data.frame to store the summary of results
# This will include metrics like total significant cis-eQTLs, unique genes, and the minimum FDR.
results_summary <- data.frame(
  strategy = character(),
  total_cis_eqtls = integer(),
  unique_genes = integer(),
  unique_snps = integer(),
  minimum_fdr = numeric(),
  stringsAsFactors = FALSE
)

# Function to execute the eQTL analysis for a given configuration
# Parameters:
#   config: A list containing the analysis name, gene expression file, and covariates file.
run_eqtl_analysis <- function(config) {
  # Define the model type for eQTL analysis
  eqtl_model <- modelLINEAR
  
  # File paths for covariates and temporary outputs
  covariates_path <- paste0(base_directory, config$covariates_file)
  temp_output_cis <- tempfile()
  temp_output_trans <- tempfile()
  
  # Print the current strategy name to track progress
  print(config$name)
  
  # Load genotype and gene expression files
  snp_file_path <- paste0(base_directory, "SNP_DEG.txt")  # SNP genotype data
  snp_location_path <- paste0(base_directory, "snpsloc_DEG.txt")  # SNP location data
  gene_expression_path <- paste0(base_directory, config$gene_expression_file)  # Gene expression data
  gene_location_path <- paste0(base_directory, "geneloc.txt")  # Gene location data
  
  # Configure SNP data
  snp_data <- SlicedData$new()
  snp_data$fileDelimiter <- "\t"  # Tab-separated file
  snp_data$fileOmitCharacters <- "NA"  # Missing data representation
  snp_data$fileSkipRows <- 1  # Skip header row
  snp_data$fileSkipColumns <- 1  # Skip ID column
  snp_data$fileSliceSize <- 20000  # Process 20,000 SNPs per slice
  snp_data$LoadFile(snp_file_path)
  
  # Configure gene expression data
  gene_data <- SlicedData$new()
  gene_data$fileDelimiter <- "\t"
  gene_data$fileOmitCharacters <- "NA"
  gene_data$fileSkipRows <- 1
  gene_data$fileSkipColumns <- 1
  gene_data$fileSliceSize <- 2000  # Process 2,000 genes per slice
  gene_data$LoadFile(gene_expression_path)
  
  # Load covariates file
  covariates_data <- SlicedData$new()
  covariates_data$fileDelimiter <- "\t"
  covariates_data$fileOmitCharacters <- "NA"
  covariates_data$fileSkipRows <- 1
  covariates_data$fileSkipColumns <- 1
  covariates_data$LoadFile(covariates_path)
  
  # Load SNP and gene position data
  snp_positions <- read.table(snp_location_path, header = TRUE, stringsAsFactors = FALSE)
  gene_positions <- read.table(gene_location_path, header = TRUE, stringsAsFactors = FALSE)
  
  # Run the Matrix_eQTL_main function to identify cis-eQTLs
  eqtl_results <- Matrix_eQTL_main(
    snps = snp_data,
    gene = gene_data,
    cvrt = covariates_data,
    output_file_name = temp_output_trans,  # Trans results output
    pvOutputThreshold = 0,  # No threshold for trans-eQTLs
    useModel = eqtl_model,
    errorCovariance = numeric(),
    verbose = TRUE,  # Print progress messages
    output_file_name.cis = temp_output_cis,  # Cis results output
    pvOutputThreshold.cis = 5e-2,  # Cis-eQTL p-value threshold
    snpspos = snp_positions,  # SNP positions
    genepos = gene_positions,  # Gene positions
    cisDist = 1e6,  # 1 Mb cis-distance
    pvalue.hist = "qqplot",  # Generate a QQ-plot
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
  )
  
  # Filter significant cis-eQTLs based on FDR < 0.05
  significant_cis_eqtls <- subset(eqtl_results$cis$eqtls, FDR < 0.05)
  
  # Calculate aggregated metrics
  total_cis_eqtls <- nrow(significant_cis_eqtls)  # Total number of significant cis-eQTLs
  unique_genes_count <- length(unique(significant_cis_eqtls$gene))  # Number of unique genes
  unique_snps_count <- length(unique(significant_cis_eqtls$snps))  # Number of unique SNPs
  minimum_fdr_value <- min(significant_cis_eqtls$FDR, na.rm = TRUE)  # Minimum FDR value
  
  # Create a summary data frame for the current strategy
  strategy_summary <- data.frame(
    strategy = config$name,
    total_cis_eqtls = total_cis_eqtls,
    unique_genes = unique_genes_count,
    unique_snps = unique_snps_count,
    minimum_fdr = minimum_fdr_value,
    stringsAsFactors = FALSE
  )
  
  # Append the summary to the global results data frame
  assign("results_summary", rbind(get("results_summary", envir = .GlobalEnv), strategy_summary), envir = .GlobalEnv)
  
  # Save results
  write.csv(significant_cis_eqtls, paste0("cis_FDR_", config$name, ".csv"))
  
  # Clean up temporary files
  unlink(temp_output_trans)
  unlink(temp_output_cis)
}

# Run the analysis for STRATEGY_A only
run_eqtl_analysis(analysis_config[[1]])

# Print the final summary table
# This table includes key metrics like total cis-eQTLs, unique genes, and minimum FDR values.
print(results_summary)


