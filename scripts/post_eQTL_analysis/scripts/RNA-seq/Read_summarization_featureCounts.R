# Read summarization with featureCounts(RSubread)
# In R
library(Rsubread)
indir <- "/mnt/sdb/projects/Iraide/MS_Iraide_Alloza_2023_10/RNA_seq/bams"
files <- list.files(path = indir, pattern = ".bam$")
empty_vec_1 <- rep("", times = length(files))
 for (i in 1:length(files)){
   file <- files[i]
   empty_vec_1[i] <- as.character(paste0(indir, "/",file))
}

indir2 <- "/mnt/sdb/projects/Iraide/MS_Iraide_Alloza_2023_10/RNA_seq/bams4"
files <- list.files(path = indir2, pattern = ".bam$")
empty_vec_2 <- rep("", times = length(files))
 for (i in 1:length(files)){
   file <- files[i]
   empty_vec_2[i] <- as.character(paste0(indir2, "/",file))
}
empty_vec_def <- c( empty_vec_1,  empty_vec_2)




# 
anno_genes <- flattenGTF("/mnt/sdb/references/Genomes/Hg38/Gencode/gencode.v44.primary_assembly.annotation.gtf.gz")
Counts_unique_reads_genes <- featureCounts(empty_vec_def, isPairedEnd = TRUE, annot.ext=anno_genes , useMetaFeatures=TRUE)
saveRDS(Counts_unique_reads_genes, file="/mnt/sdb/projects/Iraide/MS_Iraide_Alloza_2023_10/RNA_seq/counts3/Counts_unique_reads_genes_def.rds")
