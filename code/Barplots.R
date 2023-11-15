
remove(list = ls())

library(arrow)
library(readr)
library(ggplot2)

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

for (name in colnames(reduced_model)){
  print(reduced_model[[name]] == "S")
  break
}

sub_level <- colSums(reduced_model == "S")
res_level <- colSums(reduced_model == "R" | reduced_model == "I")

valid_test <- sub_level + res_level
valid_test <- valid_test[valid_test >= 50]
print(valid_test)

sub_level <- sub_level[names(valid_test)]
res_level <- res_level[names(valid_test)]

sub_prob <- (sub_level * 100) / valid_test
res_prob <- (res_level * 100) / valid_test

#sort from most to least
res_prob <- res_prob[order(res_prob)]

res_prob <- as.data.frame(res_prob)
colnames(res_prob) <- 'resisted'

#another way of sorting for data frames
#res_prob <- arrange(res_prob, resisted)

bacteria <- reorder(rownames(res_prob), res_prob$resisted)

res_ascending <- ggplot(res_prob, aes(x=reorder(rownames(res_prob), resisted), 
                                      y=resisted)) + geom_bar(stat = "identity")

res_ascending + ggtitle("Resistance Percentage by Antibiotic") +
  xlab("Antibiotic Type") + ylab("Bacteria Resisted %")

res_prob <- arrange(res_prob, desc(resisted))

res_descending <- ggplot(res_prob, aes(x=reorder(rownames(res_prob), -resisted), 
                                       y=resisted)) + geom_bar(stat = "identity")

res_descending + ggtitle("Resistance Percentage by Antibiotic") +
  xlab("Antibiotic Type") + ylab("Bacteria Resisted %")
