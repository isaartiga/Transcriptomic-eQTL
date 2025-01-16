indmiss <- read.table(file="miss.imiss", header=TRUE)
snpmiss <- read.table(file="miss.lmiss", header=TRUE)
# read data into R 

pdf("ind_miss.pdf") #indicates pdf format and gives title to file
hist(indmiss[,6], main = "Histogram individual missingness") #selects column 6, names header of file

pdf("snp_miss.pdf") 
hist(snpmiss[,5], main = "Histogram SNP missingness")  
dev.off() # shuts down the current device
