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

library(nnet)   # package needed for multinomial logistic regression
library(arrow)  # package needed to work with Parquet files
library(dplyr)
library(tidyr)
library(caret)
library(ggplot2)



# subject to change based off final implementation
amr_data_file <- "C:\\Users\\Justin\\Desktop\\School\\CIS4020\\R Scripts\\AMR_Data\\AMR Data.parquet"
drug_tiers_data_file <-"C:\\Users\\Justin\\Desktop\\School\\CIS4020\\R Scripts\\AMR_Data\\Drug Tiers.csv"


# Read the Parquet file
parquet_file <- arrow::read_parquet(amr_data_file)


# Define the predictor variables (last 57 columns)
predictor_variables <- colnames(parquet_file)[(ncol(parquet_file) - 56):ncol(parquet_file)]
predictor_variables <- c("age_year", predictor_variables)
#print(predictor_variables)


# Create a data frame for modeling
model_data <- as.data.frame(parquet_file)



reduced_model <- model_data[, predictor_variables]
# Find columns with all "NA" values
na_only_columns <- colSums(is.na(reduced_model)) == nrow(reduced_model)
#print(na_only_columns)
# Remove those columns
reduced_model <- reduced_model[, !na_only_columns, drop = FALSE]
reduced_model <- reduced_model[!is.na(reduced_model$age_year), ]
reduced_model[is.na(reduced_model)] <- "NA"

for (col in names(reduced_model)) {
  if (col != "age_year"){
    reduced_model[[col]] <- factor(reduced_model[[col]], levels = c("I", "R", "S", "N/I", "NA"))
    #print(levels(reduced_model[[col]]))
  } else {
    for (i in 1:(length(reduced_model[[col]]))){
      if (reduced_model[[col]][i] == "< 1 year"){
        reduced_model[[col]][i] <- as.integer(0.5)
      } else {
        reduced_model[[col]][i] = as.integer(reduced_model[[col]][i])
      }
    }
    
    #print(levels(reduced_model[[col]]))
  }
}
print(reduced_model$age_year)
print(levels(reduced_model$age_year))


# Check and filter columns with less than two levels
valid_columns <- character(0)
for (col in names(reduced_model)) {
  if (col == "age_year" || length(levels(reduced_model[[col]]) >= 2)) {
    valid_columns <- c(valid_columns, col)
  } else {
    cat("Column with insufficient levels:", col, "\n")
    cat("Levels:", levels(reduced_model[[col]]), "\n\n")
  }
}

# Create a subset of the dataset with valid columns
#reduced_model <- reduced_model[, c("age_year", valid_columns)]
reduced_model <- na.omit(reduced_model)
#print(reduced_model)

# Fit the logistic regression model
model <- lm(age_year ~ ., data = reduced_model)
print(summary(model))

observed_values <- reduced_model$age_year
predicted_values <- predict(model)

#print(max(na.omit(observed_values))+1)

# Create a sequence for the x-axis (e.g., 1, 2, 3, ...)
x <- seq_along(observed_values)

# Calculate the y-axis limit based on the range of numeric values in observed_values
y_max <- as.integer(max(na.omit(observed_values))) + 1

# Create a line graph
plot(x, observed_values, type = "o", col = "blue", pch = 16, ylim = c(0, y_max),
     xlab = "Test", ylab = "Age (years)", main = "Line Graph with Linear Regression Line")

# Add the regression line
lines(x, predicted_values, col = "red")
legend("topleft", legend = c("Observed", "Regression Line"), col = c("blue", "red"), lty = 1)








