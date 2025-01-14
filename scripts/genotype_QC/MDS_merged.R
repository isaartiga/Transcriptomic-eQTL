library(readr)
library(gdata)
library(dplyr)
library(data.table)

# This script generates MDS plots for population stratification.
# It includes manual adjustment for the boundary lines (abline) to ensure accurate stratification.

# Load MDS and population data
data <- read.table(file="MDS_merge2.mds", header = TRUE)
race <- read.table(file="racefile.txt", header = TRUE)
datafile <- merge(data, race, by = c("IID", "FID"))
head(datafile)

# Generate PDF plot
pdf("MDS.pdf", width = 7, height = 7)
for (i in 1:nrow(datafile)) {
    if (datafile[i, 14] == "EUR") {
        plot(datafile[i, 4], datafile[i, 5], type = "p", xlim = c(-0.1, 0.2), ylim = c(-0.15, 0.1), 
             xlab = "MDS Component 1", ylab = "MDS Component 2", pch = 1, cex = 0.5, col = "green")
    }
    par(new = TRUE)
    if (datafile[i, 14] == "ASN") {
        plot(datafile[i, 4], datafile[i, 5], type = "p", xlim = c(-0.1, 0.2), ylim = c(-0.15, 0.1), 
             xlab = "MDS Component 1", ylab = "MDS Component 2", pch = 1, cex = 0.5, col = "red")
    }
    par(new = TRUE)
    if (datafile[i, 14] == "AMR") {
        plot(datafile[i, 4], datafile[i, 5], type = "p", xlim = c(-0.1, 0.2), ylim = c(-0.15, 0.1), 
             xlab = "MDS Component 1", ylab = "MDS Component 2", pch = 1, cex = 0.5, col = "orange")
    }
    par(new = TRUE)
    if (datafile[i, 14] == "AFR") {
        plot(datafile[i, 4], datafile[i, 5], type = "p", xlim = c(-0.1, 0.2), ylim = c(-0.15, 0.1), 
             xlab = "MDS Component 1", ylab = "MDS Component 2", pch = 1, cex = 0.5, col = "blue")
    }
    par(new = TRUE)
    if (datafile[i, 14] == "OWN") {
        plot(datafile[i, 4], datafile[i, 5], type = "p", xlim = c(-0.1, 0.2), ylim = c(-0.15, 0.1), 
             xlab = "MDS Component 1", ylab = "MDS Component 2", pch = 3, cex = 0.7, col = "black")
    }
    par(new = TRUE)
}

# Add manually adjusted boundary lines
abline(v = -0.02, lty = 3) # Adjust this value manually as needed
abline(h = 0.04, lty = 3) # Adjust this value manually as needed
legend("topright", pch = c(1, 1, 1, 1, 3), c("EUR", "ASN", "AMR", "AFR", "OWN"), 
       col = c("green", "red", "orange", "blue", "black"), bty = "o", cex = 1)

dev.off()

# Generate TIFF plot
tiff("MDS.tiff", width = 400, height = 400, units = "mm", compression = "lzw", pointsize = 30, res = 400)
for (i in 1:nrow(datafile)) {
    if (datafile[i, 14] == "EUR") {
        plot(datafile[i, 4], datafile[i, 5], type = "p", xlim = c(-0.1, 0.2), ylim = c(-0.15, 0.1), 
             xlab = "MDS Component 1", ylab = "MDS Component 2", pch = 1, cex = 0.5, col = "green")
    }
    par(new = TRUE)
    if (datafile[i, 14] == "ASN") {
        plot(datafile[i, 4], datafile[i, 5], type = "p", xlim = c(-0.1, 0.2), ylim = c(-0.15, 0.1), 
             xlab = "MDS Component 1", ylab = "MDS Component 2", pch = 1, cex = 0.5, col = "red")
    }
    par(new = TRUE)
    if (datafile[i, 14] == "AMR") {
        plot(datafile[i, 4], datafile[i, 5], type = "p", xlim = c(-0.1, 0.2), ylim = c(-0.15, 0.1), 
             xlab = "MDS Component 1", ylab = "MDS Component 2", pch = 1, cex = 0.5, col = "orange")
    }
    par(new = TRUE)
    if (datafile[i, 14] == "AFR") {
        plot(datafile[i, 4], datafile[i, 5], type = "p", xlim = c(-0.1, 0.2), ylim = c(-0.15, 0.1), 
             xlab = "MDS Component 1", ylab = "MDS Component 2", pch = 1, cex = 0.5, col = "blue")
    }
    par(new = TRUE)
    if (datafile[i, 14] == "OWN") {
        plot(datafile[i, 4], datafile[i, 5], type = "p", xlim = c(-0.1, 0.2), ylim = c(-0.15, 0.1), 
             xlab = "MDS Component 1", ylab = "MDS Component 2", pch = 3, cex = 0.7, col = "black")
    }
    par(new = TRUE)
}

# Add manually adjusted boundary lines
abline(v = -0.02, lty = 3) # Adjust this value manually as needed
abline(h = 0.04, lty = 3) # Adjust this value manually as needed
legend("topright", pch = c(1, 1, 1, 1, 3), c("EUR", "ASN", "AMR", "AFR", "OWN"), 
       col = c("green", "red", "orange", "blue", "black"), bty = "o", cex = 1)

dev.off()