library(readr)
library(gdata)
library(dplyr)
library(data.table)

# This script identifies and processes cryptic relatedness between samples in GWAS data.
# It is designed to manually handle decisions regarding sample removal to optimize data retention.

# Load relatedness and missingness data
pihat_miss <- read.table("pihat.missingness.imiss", header = TRUE)
pihat <- read.table("pihat_min0.2.genome", header = TRUE)

# The following section is commented out because KING is not used here
# Uncomment if KING analysis is available
# king <- fread("king.kin0")
# merge_relatedness <- merge(pihat, king, by = c("FID1", "FID2"))
# no_merge <- anti_join(pihat, king, by = c("FID1", "FID2"))

# Add missingness data for IID1
miss <- left_join(x = pihat, y = pihat_miss[, c(2, 6)], by = c("IID1" = "IID"))
miss <- rename.vars(miss, from = "F_MISS", to = "F_MISS_1", info = TRUE)

# Add missingness data for IID2
miss <- left_join(x = miss, y = pihat_miss[, c(2, 6)], by = c("IID2" = "IID"))
miss <- rename.vars(miss, from = "F_MISS", to = "F_MISS_2", info = TRUE)

# Manually select which sample to keep based on missingness
# If F_MISS_1 > F_MISS_2, IID2 is kept; otherwise, IID1 is kept
# This step ensures manual review and input for final decisions
miss$IID <- ifelse(miss$F_MISS_1 > miss$F_MISS_2, miss$IID1, miss$IID2)
miss$FID <- ifelse(miss$F_MISS_1 > miss$F_MISS_2, miss$FID1, miss$FID2)

# Manual selection: define samples to keep (adjust these manually for each analysis)
# Replace the values in FID and IID with the manually reviewed list of samples to retain
FID <- c("KV24-1410_(Axiom_HGCoV2_1)_G03.CEL", "KV22_813.CEL", "KV24-1403_(Axiom_HGCoV2_1)_B06.CEL")
IID <- FID
miss_select <- cbind("FID" = FID, "IID" = IID)

# Prepare data for output
miss <- dplyr::select(miss, IID1, IID2, Z0, Z1, Z2, PI_HAT)

# Write outputs for manual review and subsequent use
write.table(miss, file = "pihat_info.drop", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(miss_select, file = "pihat.drop", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Notes:
# - Ensure manual input for 'FID' and 'IID' is updated before running this script.
# - This script avoids fully automated removal of samples to retain flexibility and accuracy in data selection.
