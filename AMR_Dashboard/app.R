#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(arrow)   # package needed to work with Parquet files
library(ggplot2) # package for line graph
library(nnet)   
library(dplyr)
library(tidyr)
library(caret)
library(broom)

remove(list = ls())


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("AMR Dashbaord"),

    # Sidebar with a slider input for age 
    sidebarLayout(
        sidebarPanel(
            sliderInput("given_age",
                        "What is the age of your companion patient?",
                        min = 1,
                        max = 9,
                        value = 6),
            htmlOutput("ageTextContainer")
        ),

        # Show a plot of the generated distribution
        mainPanel()
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
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
  
  # filter out the N/I as well
  tidy_data <- tidy_data %>%
    filter(!grepl("N/I", term))
  
  # only grab date of the asked for age
  filter_data <- reactive({
    filtered_data <- tidy_data[
      input$given_age >= tidy_data$estimate - tidy_data$std.error &
        input$given_age <= tidy_data$estimate + tidy_data$std.error, ]
    return(filtered_data)
  })
  
  ###############################
  # DISPLAY AGE ################
  output$ageText <- renderText({
    
    # get all matched rows
    matched_rows = filter_data()
    
    # only print it if something was found
    if (nrow(matched_rows) > 0) {
      # get the lowest p-value scored
      min_p_value_index <- which.min(matched_rows$p.value)
      matched_row <- matched_rows[min_p_value_index, , drop = FALSE]
      
      # Split the "term" value at index 2
      term_parts <- unlist(strsplit(as.character(matched_row$term), ""))
      # Construct the message
      antibiotic_result <- ""
      if (term_parts[3] == "R") {
        antibiotic_result <- "resistant"
      } else if (term_parts[3] == "S"){
        antibiotic_result <- "susceptible"
      } else {
        antibiotic_result <- "UPDATE"
      }
      output_text <- paste0(term_parts[1], term_parts[2], " has been found to be a ", antibiotic_result," antibiotic for companions of age ",input$given_age, sep = "\n")
      return(HTML(paste(output_text, collapse = "<br>"))) 
    }
    return(HTML(paste("No useful found for this age!", collapse = "<br>")))
  })
  
  output$ageTextContainer <- renderUI({
    tags$div(
      style = "padding: 10px; border: 1px solid #ccc;",
      textOutput("ageText")
    )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
