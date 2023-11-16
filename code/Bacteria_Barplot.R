
remove(list = ls())
gc()

if (length(dev.list()) != 0){
  dev.off()
}

library(arrow)
library(readr)
library(dplyr)
library(ggplot2)

amr_data_file <- "./data/AMR Data.parquet"
drug_tiers_data_file <- "./data/Drug Tiers.csv"

drug_tiers <- read_csv(drug_tiers_data_file)

parquet_file <- arrow::read_parquet(amr_data_file)
#print(unique(parquet_file$R1))

# Define the predictor variables (last 57 columns)
antibiotics <- colnames(parquet_file)[(ncol(parquet_file) - 56):ncol(parquet_file)]
bac_anti <- c("org_standard", "species", antibiotics)

bac_anti <- parquet_file[, bac_anti]
# Find columns with all "NA" values
bac_anti[bac_anti == "TF" | bac_anti == "N/I"] <- NA
na_only_columns <- colSums(is.na(bac_anti)) == nrow(bac_anti)

bac_anti <- bac_anti[, !na_only_columns, drop = FALSE]
bac_anti <- bac_anti[!is.na(bac_anti$org_standard), ]
bac_anti[is.na(bac_anti)] <- "NA"

data_rows <- length(unique(bac_anti$org_standard))
row_names <- unique(bac_anti$org_standard)
data_columns <- ncol(bac_anti[,!names(bac_anti) %in% c("org_standard", "species")])
col_names <- colnames(bac_anti[,!names(bac_anti) %in% c("org_standard", "species")])

sub_collection <- data.frame(matrix(0, data_rows, data_columns))
colnames(sub_collection) <- col_names
sub_collection$org_standard <- row_names

res_collection <- data.frame(matrix(0, data_rows, data_columns))
colnames(res_collection) <- col_names
res_collection$org_standard <- row_names

total_collection <- data.frame(matrix(0, data_rows, data_columns))
colnames(total_collection) <- col_names
total_collection$org_standard <- row_names

dog_sub <- data.frame(matrix(0, data_rows, data_columns))
colnames(dog_sub) <- col_names
dog_sub$org_standard <- row_names

dog_res <- data.frame(matrix(0, data_rows, data_columns))
colnames(dog_res) <- col_names
dog_res$org_standard <- row_names

for (i in colnames(bac_anti)){
  if (i == "org_standard" | i == "species") next
  
  anti <- bac_anti[i]
  anti$dogs <- bac_anti$species
  anti[anti == "S"] <- 1
  anti$dogs <- anti[,i] == 1 & anti$dogs == "CANINE"
  anti[anti == "R" | anti == "I" | anti == FALSE] <- 0
  anti$org_standard <- bac_anti$org_standard
  anti <- anti[anti[,i] == 1 | anti[,i] == 0,]
  sub_sums <- rowsum(as.numeric(anti[,i]), group = anti$org_standard)
  colnames(sub_sums) <- i
  dog_sums <- rowsum(as.numeric(anti[,"dogs"]), group = anti$org_standard)
  colnames(dog_sums) <- i
  
  for (j in rownames(dog_sums)){
    dog_sub[dog_sub$org_standard == j, 
                   names(dog_sub) == colnames(dog_sums)] <- (dog_sub[dog_sub$org_standard == j, 
                                                                                   names(dog_sub) == colnames(dog_sums)] + dog_sums[rownames(dog_sums) == j,])
  }
  
  for (j in rownames(sub_sums)){
    sub_collection[sub_collection$org_standard == j, 
        names(sub_collection) == colnames(sub_sums)] <- (sub_collection[sub_collection$org_standard == j, 
            names(sub_collection) == colnames(sub_sums)] + sub_sums[rownames(sub_sums) == j,])
    total_collection[total_collection$org_standard == j, 
                   names(total_collection) == colnames(sub_sums)] <- (total_collection[total_collection$org_standard == j, 
                                                                                   names(total_collection) == colnames(sub_sums)] + sub_sums[rownames(sub_sums) == j,])
  }
  
  anti <- bac_anti[i]
  anti$dogs <- bac_anti$species
  anti[anti == "R" | anti == "I"] <- 1
  anti$dogs <- anti[,i] == 1 & anti$dogs == "CANINE"
  anti[anti == "S" | anti == FALSE] <- 0
  anti$org_standard <- bac_anti$org_standard
  anti <- anti[anti[,i] == 1 | anti[,i] == 0,]
  res_sums <- rowsum(as.numeric(anti[,i]), group = anti$org_standard)
  colnames(res_sums) <- i
  dog_sums <- rowsum(as.numeric(anti[,"dogs"]), group = anti$org_standard)
  colnames(dog_sums) <- i
  
  for (j in rownames(dog_sums)){
    dog_res[dog_res$org_standard == j, 
            names(dog_res) == colnames(dog_sums)] <- (dog_res[dog_res$org_standard == j, 
                                                              names(dog_res) == colnames(dog_sums)] + dog_sums[rownames(dog_sums) == j,])
  }
  
  for (j in rownames(res_sums)){
    res_collection[res_collection$org_standard == j, 
                   names(res_collection) == colnames(res_sums)] <- (res_collection[res_collection$org_standard == j, 
                                                                                   names(res_collection) == colnames(res_sums)] + res_sums[rownames(res_sums) == j,])
    total_collection[total_collection$org_standard == j, 
                   names(total_collection) == colnames(res_sums)] <- (total_collection[total_collection$org_standard == j, 
                                                                                   names(total_collection) == colnames(res_sums)] + res_sums[rownames(res_sums) == j,])
  }

}

sub_collection$bac_sum <- rowSums(sub_collection[,!names(sub_collection) %in% c("org_standard", "species")])
res_collection$bac_sum <- rowSums(res_collection[,!names(res_collection) %in% c("org_standard", "species")])
total_collection$bac_sum <- rowSums(total_collection[,!names(total_collection) %in% c("org_standard", "species")])

res_collection$prob <- (res_collection$bac_sum / total_collection$bac_sum) * 100

res_ascending <- ggplot(res_collection, aes(x=reorder(org_standard, prob), 
                                      y=prob)) + geom_bar(stat = "identity")

res_ascending + ggtitle("Resistance Percentage by Bacteria") +
  xlab("Bacteria Type") + ylab("Antibiotics Resisted %") + ylim(0, 100)

res_descending <- ggplot(res_collection, aes(x=reorder(org_standard, -prob), 
                                            y=prob)) + geom_bar(stat = "identity")

res_descending + ggtitle("Resistance Percentage by Bacteria") +
  xlab("Bacteria Type") + ylab("Antibiotics Resisted %") + ylim(0, 100)



sub_collection$prob <- (sub_collection$bac_sum / total_collection$bac_sum) * 100

sub_ascending <- ggplot(sub_collection, aes(x=reorder(org_standard, prob), 
                                            y=prob)) + geom_bar(stat = "identity")

sub_ascending + ggtitle("Subseceptibility Percentage by Bacteria") +
  xlab("Bacteria Type") + ylab("Antibiotics Resisted %") + ylim(0, 100)

sub_descending <- ggplot(sub_collection, aes(x=reorder(org_standard, -prob), 
                                             y=prob)) + geom_bar(stat = "identity")

sub_descending + ggtitle("Subseceptibility Percentage by Bacteria") +
  xlab("Bacteria Type") + ylab("Antibiotics Resisted %") + ylim(0, 100)

#Naive Bayes Probabilities
dog_sub <- dog_sub[,!names(dog_sub) %in% c("org_standard")]
cat_sub <- sub_collection[, !names(sub_collection) %in% c("org_standard", "prob", "bac_sum")] - dog_sub
dog_res <- dog_res[,!names(dog_res) %in% c("org_standard")]
cat_res <- res_collection[, !names(res_collection) %in% c("org_standard", "prob", "bac_sum")] - dog_res

total_cats <- cat_sub + cat_res
total_dogs <- dog_sub + dog_res

total_cats[total_cats < 10] <- 0
total_dogs[total_dogs < 10] <- 0

cat_probs <- (cat_sub / total_cats) * 100
dog_probs <- (dog_sub / total_dogs) * 100

cat_probs[is.na(cat_probs)] <- 0
dog_probs[is.na(dog_probs)] <- 0

cat_probs[cat_probs == Inf] <- 0
dog_probs[dog_probs == Inf] <- 0

rownames(cat_probs) <- total_collection$org_standard
rownames(dog_probs) <- total_collection$org_standard

best_antibiotics <- function(bacteria, probs){
  if (sum(rownames(probs) == bacteria) == 0) {
    print(sprintf("Bacteria '%s' is not in this list.", bacteria))
  }
  best_options <- apply(probs[bacteria,],1,function(x) which(x==max(x)))
  best_options <- colnames(probs)[best_options]
  
  count <- 0
  for (i in best_options){
    if (count >= 3) break
    if (probs[bacteria, i] == 0) break
    print(sprintf("Antibiotic #%d: %s, Recorded Probability: %f%%", count + 1, i, probs[bacteria, i]))
    drugs <- i == drug_tiers
    for (j in 1:nrow(drugs)){
      if (drugs[j,1]){
        print(sprintf("Drug Tier: %d", as.numeric(drug_tiers[j,2])))
      }
    }
    count <- count + 1
  }
  if (count == 0){
    print(sprintf("No recorded information for bacteria '%s' on this species.", bacteria))
  }
}
print(rownames(dog_probs))
best_antibiotics("PSEUDOMONAS LUTEOLA", dog_probs)

      