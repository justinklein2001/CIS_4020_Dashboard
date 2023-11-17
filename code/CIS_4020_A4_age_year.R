###################
#
# Justin Klein
# CIS4020
# Nov 2023
#
###################
# Clear Memory
# & Load Packages
###################
remove(list=ls())
#install.packages("nnet")
#install.packages("arrow")
#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("caret")
#install.packages("ggplot2")
#install.packages("broom")

library(nnet)   
library(arrow)  
library(dplyr)
library(tidyr)
library(caret)
library(ggplot2)
library(broom)

amr_data_file <- "./data/AMR Data.parquet"

# Read in the AMR data file into memory
amr_data_file_mem <- arrow::read_parquet(amr_data_file)

# Get all the antibiotics)
required_variables <- colnames(amr_data_file_mem)[(ncol(amr_data_file_mem) - 56):ncol(amr_data_file_mem)]

# append the variable we are predicting to the data set
required_variables <- c("age_year", required_variables)

# Create a data frame to be used for modeling
model_data <- as.data.frame(amr_data_file_mem)

# Only grab the required variables from the dataset
reduced_model <- model_data[, required_variables]

# Find all columns with NA values, need to remove them
na_only_columns <- colSums(is.na(reduced_model)) == nrow(reduced_model)

# Remove the NA columns
reduced_model <- reduced_model[, !na_only_columns, drop = FALSE]

# Remove any entry that has an NA age_year entry, useless for this model
reduced_model <- reduced_model[!is.na(reduced_model$age_year), ]

# hacky solution, replace any leftover NA entries with a string so that the
# model can be computed, will be filtered out later
reduced_model[is.na(reduced_model)] <- "NA"

# start looping though the variables for pre-processing
for (col in names(reduced_model)) {
  
  # antibiotic columns
  if (col != "age_year"){
    
    # convert Antibioitc to factor to work with linear regression
    reduced_model[[col]] <- factor(reduced_model[[col]], levels = c("I", "R", "S", "N/I", "NA"))
  
  # age year col
  } else {
    
    # go through all the age_year entries
    for (i in 1:(length(reduced_model[[col]]))){
      
      # get rid of the string representation of less than a year
      if (reduced_model[[col]][i] == "< 1 year"){
        
        reduced_model[[col]][i] <- as.integer(0.5)
      
      } else {
        
        reduced_model[[col]][i] = as.integer(reduced_model[[col]][i])
      
      }
    }
  }
}

# final check for na before building model (overkill ik)
reduced_model <- na.omit(reduced_model)


# The linear regression model
model <- lm(age_year ~ ., data = reduced_model)

# Tidy up the model coefficients
tidy_data <- tidy(model)

# Filter out useless coefficients
tidy_data <- filter(tidy_data, p.value < 0.05, estimate > 0)

# this is where we actually filter out those hardcoded NA entries
tidy_data <- tidy_data %>%
  filter(!grepl("NA", term))

print(tidy_data)
# Create a plot using ggplot2
ggplot(tidy_data, aes(x = term, y = estimate)) +
  geom_point(aes(color = p.value < 0.05)) +  # Highlight significant coefficients
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0.2) +  # Confidence intervals
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") + # set the reference line to 0
  labs(title = "Linear Regression Coefficients",
       x = "Antibiotic Susceptibility",
       y = "Predicted Age") +
  theme_minimal()








