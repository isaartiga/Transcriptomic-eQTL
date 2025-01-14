# ================================
# Script: SNP to Gene Enrichment Analysis
# ================================

# ---- 1. Install and Load Required Packages ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_packages <- c("biomaRt", "dplyr", "clusterProfiler", "org.Hs.eg.db", 
                       "enrichplot", "ggplot2", "paletteer", "forcats", "stringr")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("biomaRt", "clusterProfiler", "org.Hs.eg.db", "enrichplot")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

# Load libraries
library(biomaRt)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(paletteer)
library(forcats)
library(stringr)

# ---- 2. Retrieve Genes Associated with SNPs ----
ensembl <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

# Example SNP list (replace with your real SNP list)
snp_list <- cis_FDR_A_anotado$rs  # Replace with real SNPs

results_genes <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", 
                                      "chrom_end", "associated_gene"),
                       filters = "snp_filter",
                       values = snp_list,
                       mart = ensembl)

results_genes <- results_genes %>%
  filter(chr_name %in% c(1:22)) %>%
  filter(associated_gene != "")

list_genes <- results_genes %>%
  pull(associated_gene) %>%
  strsplit(split = ",") %>%
  unlist() %>%
  unique()

# Filter Valid Gene Names 
genes_associated_valid <- bitr(list_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

genes_associated <- genes_associated_valid$SYMBOL

# ---- 3. KEGG Pathway Enrichment Analysis ----

kegg_results <- enrichKEGG(gene = genes_associated_valid$ENTREZID, 
                           organism = "hsa", 
                           pvalueCutoff = 0.05)

kegg_plot_data <- kegg_results@result %>%
  head(20) %>%
  mutate(GeneRatio = Count / as.numeric(sub("/.*", "", BgRatio)),
         Description = fct_reorder(Description, GeneRatio)) %>%
  filter(Count %in% c(6, 7))  # Filter for Count = 6 and 7

# KEGG Plot with 'Plasma' palette and consistent legends
ggplot(kegg_plot_data, aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = factor(Count), color = p.adjust)) +  # Convert Count to factor
  scale_size_manual(values = c("6" = 4, "7" = 6), name = "Counts") +  # First scale (Counts on top)
  scale_color_gradientn(colors = paletteer_c("grDevices::Plasma", 30), name = "P.adjust") +  # Second scale
  theme_minimal() +
  labs(x = "Gene Ratio", y = "", title = "Top 20 KEGG Enriched Pathways") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right")


# ---- 4. GO BP Enrichment Analysis ----
go_results <- enrichGO(gene = genes_associated, 
                       OrgDb = org.Hs.eg.db, 
                       keyType = "SYMBOL", 
                       ont = "BP", 
                       pvalueCutoff = 0.05)

go_results_simplified <- clusterProfiler::simplify(go_results, cutoff = 0.7, 
                                                   by = "p.adjust", select_fun = min)

parse_ratio <- function(ratio) {
  sapply(strsplit(ratio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
}

go_plot_data <- go_results_simplified@result %>%
  head(20) %>%
  mutate(GeneRatio = parse_ratio(GeneRatio),
         Description = fct_reorder(Description, GeneRatio))

# GO BP Plot with 'Viridis' palette and consistent legends
ggplot(go_plot_data, aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +  # Use Count as a continuous variable
  scale_size_area(max_size = 8, name = "Counts") +   # Continuous scale for sizes
  scale_color_viridis_c(option = "D", name = "P.adjust", guide = guide_colorbar(reverse = TRUE)) +  # Color scale
  theme_minimal() +
  labs(x = "Gene Ratio", y = "", title = "Top 20 GO BP Enriched Pathways") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right")


