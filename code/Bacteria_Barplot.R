
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

data_rows <- length(unique(bac_anti$org_standard))
row_names <- unique(bac_anti$org_standard)
data_columns <- ncol(bac_anti[,names(bac_anti) != "org_standard"])
col_names <- colnames(bac_anti[,names(bac_anti) != "org_standard"])

sub_collection <- data.frame(matrix(0, data_rows, data_columns))
colnames(sub_collection) <- col_names
sub_collection$org_standard <- row_names

res_collection <- data.frame(matrix(0, data_rows, data_columns))
colnames(res_collection) <- col_names
res_collection$org_standard <- row_names

total_collection <- data.frame(matrix(0, data_rows, data_columns))
colnames(total_collection) <- col_names
total_collection$org_standard <- row_names

for (i in colnames(bac_anti)){
  if (i == "org_standard") next
  
  anti <- bac_anti[i]
  anti[anti == "S"] <- 1
  anti[anti == "R" | anti == "I"] <- 0
  anti$org_standard <- bac_anti$org_standard
  anti <- anti[anti[i] == 1 | anti[i] == 0,]
  sub_sums <- rowsum(as.numeric(anti[,i]), group = anti$org_standard)
  colnames(sub_sums) <- i
  
  for (j in rownames(sub_sums)){
    sub_collection[sub_collection$org_standard == j, 
        names(sub_collection) == colnames(sub_sums)] <- (sub_collection[sub_collection$org_standard == j, 
            names(sub_collection) == colnames(sub_sums)] + sub_sums[rownames(sub_sums) == j,])
    total_collection[total_collection$org_standard == j, 
                   names(total_collection) == colnames(sub_sums)] <- (total_collection[total_collection$org_standard == j, 
                                                                                   names(total_collection) == colnames(sub_sums)] + sub_sums[rownames(sub_sums) == j,])
  }
  #print(sub_collection)
  
  anti <- bac_anti[i]
  anti[anti == "S"] <- 0
  anti[anti == "R" | anti == "I"] <- 1
  anti$org_standard <- bac_anti$org_standard
  anti <- anti[anti[i] == 1 | anti[i] == 0,]
  res_sums <- rowsum(as.numeric(anti[,i]), group = anti$org_standard)
  colnames(res_sums) <- i
  
  for (j in rownames(res_sums)){
    res_collection[res_collection$org_standard == j, 
                   names(res_collection) == colnames(res_sums)] <- (res_collection[res_collection$org_standard == j, 
                                                                                   names(res_collection) == colnames(res_sums)] + res_sums[rownames(res_sums) == j,])
    total_collection[total_collection$org_standard == j, 
                   names(total_collection) == colnames(res_sums)] <- (total_collection[total_collection$org_standard == j, 
                                                                                   names(total_collection) == colnames(res_sums)] + res_sums[rownames(res_sums) == j,])
  }
  
}

sub_collection$bac_sub <- rowSums(sub_collection[,names(sub_collection) != "org_standard"])
res_collection$bac_sub <- rowSums(res_collection[,names(res_collection) != "org_standard"])
total_collection$bac_sub <- rowSums(total_collection[,names(total_collection) != "org_standard"])

print(total_collection)

print(sub_collection$bac_sub + res_collection$bac_sub)
print(total_collection$bac_sub)


      