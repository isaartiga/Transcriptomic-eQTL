het <- read.table("het.het", header = TRUE)
pdf("heterozygosity.pdf")
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
hist(het$HET_RATE, xlab ="Heterozygosity Rate", ylab = "Frequency", main = "Heterozygosity Rate")
dev.off()
