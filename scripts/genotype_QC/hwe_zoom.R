
hwe <- read.table(file = "hwe.hwe", header = TRUE)
# read data into R 

# We grab any variants that differ in missingness between cases and controls with a p value < .00001
hwe_zoom <- hwe[hwe$P < 0.00001,]

write.table(x = hwe_zoom, file = "hwe_zoom.hwe", sep = "\t", quote = F, row.names = F)
