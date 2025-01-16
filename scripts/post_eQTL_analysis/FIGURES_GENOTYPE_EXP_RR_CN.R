# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(cowplot)  # For arranging multiple plots into grids

# 1. Reshape genotype data to long format
filtered_genotypes <- genotyped_A %>%
  filter(vcf.ID %in% test_RRvsCN_genotype_A_pval$SNP)

genotype_long <- filtered_genotypes %>%
  pivot_longer(
    cols = -vcf.ID,  # Exclude vcf.ID column
    names_to = "Sample",
    values_to = "Genotype"
  )

# 2. Join with annotated variant information
annotated_genotype <- genotype_long %>%
  inner_join(cis_FDR_A_annotated, by = c("vcf.ID" = "rs")) %>%
  mutate(
    Genotype_Label = case_when(
      Genotype == 0 ~ paste(Ref.allele, Ref.allele, sep = ""),
      Genotype == 1 ~ paste(Ref.allele, Alt.allele, sep = ""),
      Genotype == 2 ~ paste(Alt.allele, Alt.allele, sep = ""),
      TRUE ~ NA_character_
    )
  )

# 3. Reshape expression data to long format
expression_long <- expression %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "Expression"
  )

# 4. Combine genotype and expression data
merged_data <- annotated_genotype %>%
  inner_join(expression_long, by = c("Gene", "Sample")) %>%
  filter(!is.na(Genotype_Label))  # Remove rows with missing genotype labels

# Add a Group column to differentiate Control and RR samples
merged_data <- merged_data %>%
  mutate(
    Group = case_when(
      grepl("^CN", Sample) ~ "Control",
      grepl("^RR", Sample) ~ "RR",
      TRUE ~ "Other"
    )
  )

# Unique SNP-Gene pairs
unique_snp_gene_pairs <- unique(merged_data[, c("vcf.ID", "Gene")])

# Generate plots for each SNP-Gene pair
plots_RR_CN_A <- unique_snp_gene_pairs %>%
  pmap(function(vcf.ID, Gene) {
    # Filter data for the current SNP-Gene pair
    pair_data <- merged_data %>% filter(vcf.ID == !!vcf.ID, Gene == !!Gene)
    
    # Define specific genotype order
    valid_levels <- c(
      paste(unique(pair_data$Ref.allele), unique(pair_data$Ref.allele), sep = ""),
      paste(unique(pair_data$Ref.allele), unique(pair_data$Alt.allele), sep = ""),
      paste(unique(pair_data$Alt.allele), unique(pair_data$Alt.allele), sep = "")
    )
    pair_data <- pair_data %>%
      mutate(Genotype_Label = factor(Genotype_Label, levels = valid_levels))
    
    # Create the plot with specific colors and stronger points
    ggplot(pair_data, aes(x = Genotype_Label, y = Expression, color = Group)) +
      geom_jitter(width = 0.2, size = 2, alpha = 0.9) +  # Reduce transparency (stronger points)
      scale_color_manual(values = c("Control" = "#66CDAA", "RR" = "#EE9572", "Other" = "grey")) +
      theme_minimal() +
      labs(
        x = "Genotype",
        y = "Expression Level",
        title = paste("eQTL Analysis\nSNP:", vcf.ID, "Gene:", Gene),
        color = "Group"
      ) +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  })

# Split plots into batches of 16 per page
plots_per_page <- 16
plot_batches <- split(plots_RR_CN_A, ceiling(seq_along(plots_RR_CN_A) / plots_per_page))

# Generate and save each page of plots
output_dir <- "plots_RRvsCN"  # Define output directory
if (!dir.exists(output_dir)) dir.create(output_dir)  # Create directory if it doesn't exist

for (i in seq_along(plot_batches)) {
  combined_plot <- plot_grid(plotlist = plot_batches[[i]], ncol = 4, nrow = 4)
  
  # Display the current page
  print(combined_plot)
  
  # Save the page as a file
  ggsave(
    filename = paste0(output_dir, "/eQTL_A_page_", i, ".png"),
    plot = combined_plot,
    width = 16,
    height = 16
  )
}


