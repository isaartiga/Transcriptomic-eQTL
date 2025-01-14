### eQTL Analysis: Combining Genotype and Expression Data for Visualization
# Author: Isabel Artiga Folch
# Date: 17/12/2024
# Description: This script transforms genotype and expression datasets, merges them, and creates violin plots 
#              for each SNP-Gene pair to analyze expression quantitative trait loci (eQTL).

# Load necessary libraries
library(tidyverse)  # Data manipulation and visualization
library(forcats)    # Factor manipulation
library(cowplot)    # Combine multiple plots

# 1. Transform genotype data to long format
# Convert the genotype dataset from wide to long format for analysis
# 'vcf.ID' identifies SNPs, and other columns represent sample genotypes

genotipado_long <- genotipado_A %>%
  pivot_longer(
    cols = -vcf.ID,  # Exclude 'vcf.ID' (SNP identifier) from pivot
    names_to = "Sample",       # Rename columns to 'Sample'
    values_to = "Genotype"     # Place genotype values in a column named 'Genotype'
  )

# 2. Merge with annotated variant information
# Add reference and alternate allele information using the annotation dataset

genotipo_annotado <- genotipado_long %>%
  inner_join(cis_FDR_A_anotado, by = c("vcf.ID" = "rs")) %>%  # Merge on SNP ID (vcf.ID)
  mutate(
    # Create a new column for labeled genotypes
    Genotype_Label = case_when(
      Genotype == 0 ~ paste(Ref.allele, Ref.allele, sep = ""),  # Homozygous reference
      Genotype == 1 ~ paste(Ref.allele, Alt.allele, sep = ""),  # Heterozygous
      Genotype == 2 ~ paste(Alt.allele, Alt.allele, sep = ""),  # Homozygous alternate
      TRUE ~ NA_character_   # Default to NA for missing values
    )
  )

# 3. Transform expression data to long format
# Convert the gene expression dataset from wide to long format

expression_long <- expression %>%
  pivot_longer(
    cols = -Gene,             # Exclude 'Gene' (gene identifier)
    names_to = "Sample",     # Rename columns to 'Sample'
    values_to = "Expression" # Place expression values in 'Expression'
  )

# 4. Combine genotype and expression data
# Merge genotype and expression data based on Gene and Sample identifiers

merged_data <- genotipo_annotado %>%
  inner_join(expression_long, by = c("Gene", "Sample"))

# Convert to data frame for compatibility
merged_data <- as.data.frame(merged_data)

# Extract unique SNP-Gene pairs
unique_snp_gene_pairs <- unique(merged_data[, c("vcf.ID", "Gene")])

# Filter out rows with missing genotype labels
merged_data <- merged_data %>%
  filter(!is.na(Genotype_Label))  # Remove rows where Genotype_Label is NA

# 5. Create violin plots for each SNP-Gene pair
# Loop through each unique SNP-Gene combination to create visualizations

plots_A <- unique_snp_gene_pairs %>%
  pmap(function(vcf.ID, Gene) {
    # Filter data for the current SNP-Gene pair
    pair_data <- merged_data %>% filter(vcf.ID == !!vcf.ID, Gene == !!Gene)
    
    # Get the possible levels that are present in the current data
    valid_levels <- c(
      paste(unique(pair_data$Ref.allele), unique(pair_data$Ref.allele), sep = ""), 
      paste(unique(pair_data$Ref.allele), unique(pair_data$Alt.allele), sep = ""), 
      paste(unique(pair_data$Alt.allele), unique(pair_data$Alt.allele), sep = "")
    )
    
    # Reorder Genotype_Label if level matches
    pair_data <- pair_data %>%
      mutate(Genotype_Label = factor(Genotype_Label, levels = valid_levels)) %>%
      mutate(Genotype_Label = fct_infreq(Genotype_Label))  # Reorder by frequency
    
    # Create the plot excluding NA genotypes
    ggplot(pair_data, aes(x = Genotype_Label, y = Expression, fill = Genotype_Label)) +
      geom_violin(trim = TRUE, alpha = 0.5) +
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 0.5, alpha = 0.7) +
      theme_minimal() +
      labs(
        x = "Genotype",
        y = "Expression Level",
        title = paste("eQTL Analysis\nSNP:", vcf.ID, "Gene:", Gene)
      ) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
      )
  })


# 6. Split plots into pages and save them
# Divide the generated plots into batches of 16 per page

plots_per_page <- 16
plot_batches <- split(plots_A, ceiling(seq_along(plots_A) / plots_per_page))

# Loop through each batch to combine and save the plots
for (i in seq_along(plot_batches)) {
  # Combine plots into a grid
  combined_plot <- plot_grid(plotlist = plot_batches[[i]], ncol = 4, nrow = 4)
  
  # Display the current page of plots
  print(combined_plot)
  
  # Save the page as a PNG file
  ggsave(
    filename = paste0("violin_plots_2.0/eQTL_A_page_", i, ".png"),
    plot = combined_plot,
    width = 16,  # Width in inches
    height = 16  # Height in inches
  )
}

### End of Script
# Notes:
# - Input files: 'genotipado_A', 'cis_FDR_A_anotado', 'expression'
# - Output files: Multiple PNG files (one per page of plots)
# - Ensure input data are properly formatted before running the script.
