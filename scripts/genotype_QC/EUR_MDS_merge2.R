library(readr)
library(gdata)
library(dplyr)
library(data.table)

# This script generates MDS plots and allows filtering of EUR samples.
# Boundary limits can be determined manually after generating the plot.

# Load MDS and population data
data <- read.table(file = "MDS_merge2.mds", header = TRUE)
race <- read.table(file = "racefile.txt", header = TRUE)
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
legend("topright", pch = c(1, 1, 1, 1, 3), c("EUR", "ASN", "AMR", "AFR", "OWN"), 
       col = c("green", "red", "orange", "blue", "black"), bty = "o", cex = 1)

# Manual step: Review the plot (MDS_no_lines.pdf) to determine cutoff limits
# Set the cutoff limits based on the plot
cutoff_C1 <- -0.02  # Example: Adjust this after reviewing the plot
cutoff_C2 <- 0.04   # Example: Adjust this after reviewing the plot

dev.off()


# Filter EUR samples based on determined cutoffs
data_EUR <- datafile[datafile$C1 <= cutoff_C1 & datafile$C2 >= cutoff_C2, c("IID", "FID")]

# Save the filtered EUR file
write.table(data_EUR, "EUR_MDS_merge2", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

tiff("MDS.tiff", width=400,height=400,units="mm",compression=c("lzw"),pointsize=30, res=400)
for (i in 1:nrow(datafile))
{
    if (datafile[i,14]=="EUR") {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="green")}
    par(new=T)
    if (datafile[i,14]=="ASN") {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="red")}
    par(new=T)
    if (datafile[i,14]=="AMR") {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="orange")}
    par(new=T)
    if (datafile[i,14]=="AFR") {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="blue")}
    par(new=T)
    if (datafile[i,14]=="OWN") {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=3,cex=0.7,col="black")}
    par(new=T)
}

abline(v=cutoff_C1,lty=3)
abline(h=cutoff_C2,lty=3)
legend("topright", pch=c(1,1,1,1,3),c("EUR","ASN","AMR","AFR","OWN"),col=c("green","red","orange","blue","black"),bty="o",cex=1)

dev.off()
