
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
print(unique(parquet_file$R1))

# Define the predictor variables (last 57 columns)
antibiotics <- colnames(parquet_file)[(ncol(parquet_file) - 56):ncol(parquet_file)]

antibiotics <- parquet_file[, antibiotics]
# Find columns with all "NA" values
na_only_columns <- colSums(is.na(antibiotics)) == nrow(antibiotics)

antibiotics <- antibiotics[, !na_only_columns, drop = FALSE]
antibiotics[is.na(antibiotics)] <- "NA"

sub_level <- colSums(antibiotics == "S")
res_level <- colSums(antibiotics == "R" | antibiotics == "I")

valid_test <- sub_level + res_level
valid_test <- valid_test[valid_test >= 500]

sub_level <- sub_level[names(valid_test)]
res_level <- res_level[names(valid_test)]

sub_prob <- (sub_level * 100) / valid_test
res_prob <- (res_level * 100) / valid_test

#sort from most to least
#res_prob <- res_prob[order(res_prob)]

res_prob <- as.data.frame(res_prob)
colnames(res_prob) <- 'resisted'

#another way of sorting for data frames
#res_prob <- arrange(res_prob, resisted)

res_ascending <- ggplot(res_prob, aes(x=reorder(rownames(res_prob), resisted), 
                                      y=resisted)) + geom_bar(stat = "identity")

res_ascending + ggtitle("Resistance Percentage by Antibiotic") +
  xlab("Antibiotic Type") + ylab("Bacteria Resisted %") + ylim(0, 100)

res_descending <- ggplot(res_prob, aes(x=reorder(rownames(res_prob), -resisted), 
                                       y=resisted)) + geom_bar(stat = "identity")

res_descending + ggtitle("Resistance Percentage by Antibiotic") +
  xlab("Antibiotic Type") + ylab("Bacteria Resisted %") + ylim(0, 100)

#sort from most to least
#sub_prob <- sub_prob[order(sub_prob)]

sub_prob <- as.data.frame(sub_prob)
colnames(sub_prob) <- 'subsceptible'

#another way of sorting for data frames
#res_prob <- arrange(res_prob, resisted)

sub_ascending <- ggplot(sub_prob, aes(x=reorder(rownames(sub_prob), subsceptible), 
                                      y=subsceptible)) + geom_bar(stat = "identity")

sub_ascending + ggtitle("Subsceptibility Percentage by Antibiotic") +
  xlab("Antibiotic Type") + ylab("Bacteria Subsceptibility %") + ylim(0, 100)

sub_descending <- ggplot(sub_prob, aes(x=reorder(rownames(sub_prob), -subsceptible), 
                                      y=subsceptible)) + geom_bar(stat = "identity")

sub_descending + ggtitle("Subsceptibility Percentage by Antibiotic") +
  xlab("Antibiotic Type") + ylab("Bacteria Subsceptibility %") + ylim(0, 100)
