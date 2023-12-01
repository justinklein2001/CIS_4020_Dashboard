
remove(list = ls())

library(arrow)
library(readr)

amr_data_file <- "./data/AMR Data.parquet"
drug_tiers_data_file <- "./data/Drug Tiers.csv"

parquet_file <- arrow::read_parquet(amr_data_file)
View(parquet_file)
print(unique(parquet_file$R1))

# Define the predictor variables (last 57 columns)
predictor_variables <- colnames(parquet_file)[(ncol(parquet_file) - 56):ncol(parquet_file)]

reduced_model <- parquet_file[, predictor_variables]
# Find columns with all "NA" values
na_only_columns <- colSums(is.na(reduced_model)) == nrow(reduced_model)

reduced_model <- reduced_model[, !na_only_columns, drop = FALSE]
reduced_model[is.na(reduced_model)] <- "NA"

View(reduced_model)

print(colSums(reduced_model == "S"))
