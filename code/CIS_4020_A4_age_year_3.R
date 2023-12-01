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
#install.packages("tidyverse")
library(nnet)   
library(arrow)  
library(dplyr)
library(tidyr)
library(caret)
library(ggplot2)
library(tidyverse)
library(arrow)

# Load AMR data
amr_data_file <- "./data/AMR Data.parquet"
amr_data_file_mem <- arrow::read_parquet(amr_data_file)

# Get antibiotics columns
antibiotic_columns <- colnames(amr_data_file_mem)[(ncol(amr_data_file_mem) - 56):ncol(amr_data_file_mem)]

# Append the variable we are predicting to the data set
model_data <- amr_data_file_mem %>%
  select(age_year, all_of(antibiotic_columns)) %>%
  mutate(age_year = ifelse(age_year == "< 1 year", 0.5, as.numeric(age_year)),
         age_group = ifelse(age_year < 2, "young", "old")) %>%
  mutate(across(antibiotic_columns, factor, levels = c("I", "R", "S", "N/I", "NA")),
         age_group = factor(age_group, levels = c("young", "old")))  # Ensure age_group has two levels

print(model_data)

# Find all columns with NA values, need to remove them
na_only_columns <- colSums(is.na(model_data)) == nrow(model_data)

# Remove the NA columns
model_data <- model_data[, !na_only_columns, drop = FALSE]

# Remove any entry that has an NA age_year entry, useless for this model
model_data <- model_data[!is.na(model_data$age_year), ]


# hacky solution, replace any leftover NA entries with a string so that the
# model can be computed, will be filtered out later
model_data[is.na(model_data)] <- "NA"

print(model_data)
# Logistic regression model
logistic_model <- glm(age_group ~ ., data = model_data, family = "binomial")

# Print summary of the model
summary(logistic_model)

# Make predictions
predictions <- predict(logistic_model, newdata = model_data, type = "response")

# Assuming 'actual_values' is the actual values of the dependent variable
actual_values <- ifelse(model_data$age_group == "young", 1, 0)

# Calculate predicted probabilities
predicted_probabilities <- ifelse(predictions > 0.5, 1, 0)

# Print or analyze predictions
print(predicted_probabilities)

# Assess model accuracy
accuracy <- mean(predicted_probabilities == actual_values)

# Print or use accuracy metric
cat("Accuracy:", accuracy, "\n")

# Create a plot using ggplot2
tidy_data <- broom::tidy(logistic_model)
ggplot(tidy_data, aes(x = term, y = estimate)) +
  geom_point(aes(color = p.value < 0.05)) +  # Highlight significant coefficients
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0.2) +  # Confidence intervals
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") + # set the reference line to 0
  labs(title = "Logistic Regression Coefficients",
       x = "Antibiotic Susceptibility",
       y = "Predicted Age Group") +
  theme_minimal()

                       