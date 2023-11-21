#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#install.packages("shiny")
#install.packages("shinydashboard")

library(shiny)
library(shinydashboard)
remove(list = ls())

if (length(dev.list()) != 0){
  dev.off()
}

library(arrow)
library(readr)
library(dplyr)
library(ggplot2)
library(nnet)   
library(tidyr)
library(caret)
library(broom)

amr_data_file <- "./data/AMR Data.parquet"
drug_tiers_data_file <- "./data/Drug Tiers.csv"

drug_tiers <- read_csv(drug_tiers_data_file)

parquet_file <- arrow::read_parquet(amr_data_file)

####################

## JUSTIN START HERE
# store counts of each type of year in a table
test_counts <- table(parquet_file[["order_month"]])

# make counts a df
months_count_df <- data.frame(value = as.numeric(names(test_counts)), count = as.numeric(test_counts))

# Get all the antibiotics)
required_variables <- colnames(parquet_file)[(ncol(parquet_file) - 56):ncol(parquet_file)]

# append the variable we are predicting to the data set
required_variables <- c("age_year", required_variables)

# Create a data frame to be used for modeling
model_data <- as.data.frame(parquet_file)

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


## JUSTIN END HERE
#################

# Define the predictor variables (last 57 columns)
antibiotics <- colnames(parquet_file)[(ncol(parquet_file) - 56):ncol(parquet_file)]
bac_anti <- c("org_standard", "species", antibiotics)

bac_anti <- parquet_file[, bac_anti]

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
colnames(res_prob) <- 'chance'

sub_prob <- as.data.frame(sub_prob)
colnames(sub_prob) <- 'chance'

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
sub_collection$prob <- (sub_collection$bac_sum / total_collection$bac_sum) * 100

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

best_antibiotics <- function(bacteria, species){
  if (species == "Canine"){
    probs <- dog_probs
  }
  else{
    probs <- cat_probs
  }
  if (sum(rownames(probs) == bacteria) == 0) {
    print(sprintf("Bacteria '%s' is not in this list.", bacteria))
  }
  
  best_options <- t(probs[bacteria,])
  colnames(best_options) <- "bacteria"
  best_options <- best_options[order(best_options, decreasing = TRUE), ,drop = FALSE]
  
  best_labels <- rownames(best_options)
  
  best_options <- best_options[1:3]
  best_labels <- best_labels[1:3]
  
  if (best_options[1] == 0){
    return(sprintf("No recorded information for bacteria '%s' on this species.", bacteria))
  }
  text_return <- ""
  for (i in 1:3){
    text_return <- paste(text_return, sprintf("Antibiotic #%d: %s, Recorded Probability: %f%%", i, best_labels[i], best_options[i]), sep="\n")
    drugs <- best_labels[i] == drug_tiers
    for (j in 1:nrow(drugs)){
      if (drugs[j,1]){
        text_return <- paste(text_return, sprintf("Drug Tier: %d", as.numeric(drug_tiers[j,2])), sep="\n")
      }
    }
  }
  text_return
}
cat(best_antibiotics("E COLI", "Canine"))

# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "AMR Dashboard - Group Burnout"),
  dashboardSidebar(
    sidebarMenu(
      menuItem(" Home", tabName = "home", icon = icon("home")),
      menuItem(" About", tabName = "about", icon = icon("info"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "home",
              fluidRow(
                box(
                  title = "Welcome message",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 12,
                  h3("Welcome back! Please tell us about your patient:")
                ),
              ),
              fluidRow(align = "center", label="Welcome back! Tell us about your patient:",
                       box(width = 2, offset = 0, align = "center", height="auto",
                           selectInput(
                             inputId = "species",
                             label = "What species is your patient?",
                             choices = c("Canine", "Feline"),
                             selected = "Canine",
                             multiple = FALSE,
                             selectize = TRUE,
                           ),
                           textOutput("text"),
                       ),
                       box(width = 3, offset = 4, align = "center", height=100,
                           selectizeInput(
                             inputId = "search", 
                             label = "What's the observed bacteria in your patient?",
                             multiple = FALSE,
                             choices = c("Best Antibiotics Search" = "", sort(total_collection$org_standard)),
                             options = list(
                               create = FALSE,
                               placeholder = "Input Bacteria Here",
                               maxItems = '1',
                               onDropdownOpen = I("function($dropdown) {if (!this.lastQuery.length) {this.close(); this.settings.openOnFocus = false;}}"),
                               onType = I("function (str) {if (str === \"\") {this.close();}}")
                             )
                           ),
                       ),
                       box(width = 4, align = "center", height=100,
                              sliderInput("given_age",
                                          "What is the age of your companion patient?",
                                          min = 1,
                                          max = 9,
                                          value = 6),
                       ),
                       box(width = 3, align = "center", height=100,
                           htmlOutput("ageTextContainer"),
                       ),
              ),
              fluidRow(align = "center",
                       verbatimTextOutput("best_antibiotics"),
              ),
              fluidRow(
                box(
                  title = "Need more information?",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 12,
                  h4("Still haven't come to a consensus? Feel free to explore more general data to reach a diagnosis!")
                ),
              ),
              fluidRow(
                box(
                  title = "Explore Antibiotic Resistance and Susceptibility levels:",
                  status = "primary",
                  plotOutput("antiPlot"), width = 6),
                box(
                  title = "Explore Bacteria Resistance and Susceptibility levels:",
                  status = "primary",
                  plotOutput("bacPlot"), width = 6)
              ),
              fluidRow(
                box(width = 4, offset = 0, align = "center", height = 150,
                       sliderInput("bins",
                                   "# of antibiotics to display:",
                                   min = 4,
                                   max = nrow(res_prob),
                                   value = 20,
                                   step = 1,
                                   width = "75%"),
                ),
                box(width = 2, offset = 0, align = "center", height = 150, title = "Modify Graph:",
                       actionButton(
                         "swapplot_anti", 
                         "Subsceptibility",
                         style='padding:4px; font-size:80%; margin-top: 15%; width: 80px;'
                       ),
                       actionButton(
                         "sorting_anti", 
                         "Ascending",
                         style='padding:4px; font-size:80%; margin-top: 15%; width: 80px;'
                       )
                ),
                box(width = 4, offset = 0, align = "center", height = 150,
                       sliderInput("bins_bac",
                                   "# of bacteria to display:",
                                   min = 2,
                                   max = 6,#nrow(res_collection) - 1,
                                   value = 4,
                                   step = 1,
                                   width = "75%"),
                ),
                box(width = 2, offset = 0,
                       align = "center", height = 150, title = "Modify Graph:",
                       actionButton(
                         "swapplot_bac", 
                         "Subsceptibility",
                         style='padding:4px; font-size:80%; margin-top: 15%; width: 80px;'
                       ),
                       actionButton(
                         "sorting_bac", 
                         "Ascending",
                         style='padding:4px; font-size:80%; margin-top: 15%; width: 80px;'
                       )
                ),
                fluidRow(align = "center",
                         box(width = 8, align = "center", status = "primary",
                             plotOutput("num_tests"),
                         ),
                         box(width = 4, align = "center", status = "primary",
                             "Based off data from the last 4 years, it looks like November has a low amount of AMR testing.",
                             box(
                               title = "Current Saffing Needs:",
                               status = "primary",
                               solidHeader = TRUE,
                               collapsible = TRUE,
                               width = 12,
                               h4("Low amount of staff required till EOY. Give them some time off!")
                             ),
                         ),
                ),
              ),
      ),
      tabItem(tabName = "about",
              fluidRow(
                box("Content for the about tab", width = 12)
              )
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  params_anti <- reactiveValues(
    prob = arrange(res_prob, desc(chance)), 
    descending = TRUE, 
    res = TRUE
  )
  
  observeEvent(input$sorting_anti, {
    if (params_anti$descending){
      params_anti$prob <- arrange(params_anti$prob, chance)
      params_anti$descending <- FALSE
      updateActionButton(session, "sorting_anti", label = "Descending")
    }
    else{
      params_anti$prob <- arrange(params_anti$prob, desc(chance))
      params_anti$descending <- TRUE
      updateActionButton(session, "sorting_anti", label = "Ascending")
    }
  })
  
  observeEvent(input$swapplot_anti, {
    if (params_anti$res){
      params_anti$prob <- sub_prob
      params_anti$res <- FALSE
      updateActionButton(session, "swapplot_anti", label = "Subsceptibility")
    }
    else{
      params_anti$prob <- res_prob
      params_anti$res <- TRUE
      updateActionButton(session, "swapplot_anti", label = "Resistance")
    }
    
    if (!params_anti$descending){
      params_anti$prob <- arrange(params_anti$prob, chance)
    }
    else{
      params_anti$prob <- arrange(params_anti$prob, desc(chance))
    }
  })
  
  output$antiPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    bins <- as.data.frame(params_anti$prob[0:input$bins,, drop = FALSE])
    
    if (params_anti$descending){
      plot <- ggplot(bins, aes(x=reorder(rownames(bins), -chance), 
                               y=chance)) + geom_bar(stat = "identity",fill = "#3c8dbc")
    }
    else{
      plot <- ggplot(bins, aes(x=reorder(rownames(bins), chance), 
                               y=chance)) + geom_bar(stat = "identity",fill = "#3c8dbc")
    }
    
    if (params_anti$res){
      plot + ggtitle("Resistance Percentage by Antibiotic") +
        xlab("Antibiotic Type") + ylab("Bacteria Resistance %") + ylim(0, 100)
    }
    else{
      plot + ggtitle("Subsceptibility Percentage by Antibiotic") +
        xlab("Antibiotic Type") + ylab("Bacteria Subsceptibility %") + ylim(0, 100)
    }
  })
  
  params_bac <- reactiveValues(
    prob = arrange(res_collection, desc(prob)), 
    descending = TRUE, 
    res = TRUE
  )
  
  observeEvent(input$sorting_bac, {
    if (params_bac$descending){
      params_bac$prob <- arrange(params_bac$prob, prob)
      params_bac$descending <- FALSE
      updateActionButton(session, "sorting_bac", label = "Descending")
    }
    else{
      params_bac$prob <- arrange(params_bac$prob, desc(prob))
      params_bac$descending <- TRUE
      updateActionButton(session, "sorting_bac", label = "Ascending")
    }
  })
  
  observeEvent(input$swapplot_bac, {
    if (params_bac$res){
      params_bac$prob <- sub_collection
      params_bac$res <- FALSE
      updateActionButton(session, "swapplot_bac", label = "Subsceptibility")
    }
    else{
      params_bac$prob <- res_collection
      params_bac$res <- TRUE
      updateActionButton(session, "swapplot_bac", label = "Resistance")
    }
    
    if (!params_bac$descending){
      params_bac$prob <- arrange(params_bac$prob, prob)
    }
    else{
      params_bac$prob <- arrange(params_bac$prob, desc(prob))
    }
  })
  
  output$bacPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    bins_bac <- as.data.frame(params_bac$prob[0:input$bins_bac,, drop = FALSE])
    
    if (params_bac$descending){
      plot <- ggplot(bins_bac, aes(x=reorder(org_standard, -prob), 
                                   y=prob)) + geom_bar(stat = "identity", fill = "#3c8dbc")
    }
    else{
      plot <- ggplot(bins_bac, aes(x=reorder(org_standard, prob), 
                                   y=prob)) + geom_bar(stat = "identity", fill = "#3c8dbc")
    }
    
    if (params_bac$res){
      plot + ggtitle("Resistance Percentage by Bacteria") +
        xlab("Bacteria Type") + ylab("Bacteria Resistance %") + ylim(0, 100)
    }
    else{
      plot + ggtitle("Subsceptibility Percentage by Bacteria") +
        xlab("Bacteria Type") + ylab("Bacteria Subsceptibility %") + ylim(0, 100)
    }
  })
  output$best_antibiotics <- renderText({
    if (input$search != "") best_antibiotics(input$search, input$species)
  })
  
  ## JUSTIN STUFF HERE
  
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
      output_text <- paste0(term_parts[1], term_parts[2], " can be found to be a ", antibiotic_result," antibiotic for companions of age ",input$given_age, sep = "\n")
      return(HTML(paste(output_text, collapse = "<br>"))) 
    }
    return(HTML(paste("No useful information found for this age!", collapse = "<br>")))
  })
  
  output$ageTextContainer <- renderUI({
    tags$div(
      style = "padding: 10px; border: 1px solid #ccc;",
      textOutput("ageText")
    )
  })
  
  ###############################
  # DISPLAY TESTS ################
  output$num_tests <- renderPlot({
    # plot that baby
    ggplot(months_count_df, aes(x = value, y = count)) +
      geom_line(color = "#3c8dbc") +
      geom_point(color = "#3c8dbc", size = 3) +
      labs("Number of Tests Per Month",
           x = "Months",
           y = "Number of Tests",
           title = "Number of AMR Tests per Month") +
      scale_x_continuous(breaks = 1:12)
  })
  ## JUSTIN END
}

# Run the application
shinyApp(ui, server)
  
