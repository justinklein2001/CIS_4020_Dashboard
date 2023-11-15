
remove(list = ls())

if (length(dev.list()) != 0){
  dev.off()
}

library(arrow)
library(readr)
library(dplyr)
library(ggplot2)

amr_data_file <- "./data/AMR Data.parquet"
drug_tiers_data_file <- "./data/Drug Tiers.csv"

parquet_file <- arrow::read_parquet(amr_data_file)
#print(unique(parquet_file$R1))

# Define the predictor variables (last 57 columns)
antibiotics <- colnames(parquet_file)[(ncol(parquet_file) - 56):ncol(parquet_file)]
bac_anti <- c("org_standard", antibiotics)

bac_anti <- parquet_file[, bac_anti]
# Find columns with all "NA" values
na_only_columns <- colSums(is.na(bac_anti)) == nrow(bac_anti)

bac_anti <- bac_anti[, !na_only_columns, drop = FALSE]
bac_anti <- bac_anti[!is.na(bac_anti$org_standard), ]
bac_anti[is.na(bac_anti)] <- "NA"

num_bacteria = n_distinct(bac_anti$org_standard)
print(num_bacteria)

temp <- bac_anti[,!(names(bac_anti) == "org_standard")]

bac_sub_level <- rowSums(temp == "S")

bac_anti$bac_sub <- bac_sub_level
print(bac_anti)

#TO BE DONE BELOW

bac_anti <- rowsum(bac_anti, group = colnames(bac_anti), na.rm = T)

sub_level <- colSums(bac_anti == "S")
res_level <- colSums(bac_anti == "R" | bac_anti == "I")

valid_test <- sub_level + res_level
valid_test <- valid_test[valid_test >= 500]
print(valid_test)

dat <- data.frame(
  value = runif(100),
  group = sample(1:20, 100, replace = T)
)
print(dat)
rowsum(dat$value, dat$group)

      