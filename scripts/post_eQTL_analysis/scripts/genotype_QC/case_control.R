case_control <- read.table(file = "case-control.missing", header=TRUE)
# read data into R 

# We grab any variants that differ in missingness between cases and controls with a p value < .00001
case_control_drop <- case_control[case_control$P < 0.00001,]

write.table(x = case_control_drop,file = "case_control.drop", sep = "\t", quote = F, row.names = F)

