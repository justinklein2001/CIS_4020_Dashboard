#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# un-comment these lines to install packages
#install.packages("shiny")
#install.packages("shinydashboard")
#install.packages("arrow")
#install.packages("readr")
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("nnet")   
#install.packages("tidyr")
#install.packages("caret")
#install.packages("broom")

remove(list = ls())

if (length(dev.list()) != 0){
  dev.off()
}

library(shiny)
library(shinydashboard)
library(arrow)
library(readr)
library(dplyr)
library(ggplot2)
library(nnet)   
library(tidyr)
library(caret)
library(broom)


# reading in AMR data files (should be platform independent)
amr_data_file <-  file.path(getwd(), "../data/AMR Data.parquet")
drug_tiers_data_file <- file.path(getwd(), "../data/Drug Tiers.csv")

drug_tiers <- read_csv(drug_tiers_data_file)

parquet_file <- arrow::read_parquet(amr_data_file)

#########################################
## LINEAR REGRESSION PRE-PROCESSING ####
#######################################

# store counts of each type of year in a table
test_counts <- table(parquet_file[["order_month"]])

# make counts a df (NOT FOR LIN REGRESSION, but need to do this here)
months_count_df <- data.frame(value = as.numeric(names(test_counts)), count = as.numeric(test_counts))

# Get all the antibiotics
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


#########################################
## NAIVE BAYES PRE-PROCESSING ##########
#######################################

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

rownames(total_dogs) <- total_collection$org_standard
rownames(total_cats) <- total_collection$org_standard

best_antibiotics <- function(bacteria, species){
  if (species == "Canine"){
    probs <- dog_probs
    counts <- total_dogs
  }
  else{
    probs <- cat_probs
    counts <- total_cats
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
  total_pets <- counts[bacteria, best_labels]
  
  if (best_options[1] == 0){
    return(sprintf("No recorded information for bacteria '%s' on this species.", bacteria))
  }
  text_return <- ""
  for (i in 1:3){
    text_return <- paste(text_return, sprintf("Antibiotic #%d: %s, Recorded Probability: %f%%", i, best_labels[i], best_options[i]), sep="\n")
    text_return <- paste(text_return, sprintf("Total Tests Done: %d", counts[bacteria, best_labels[i]]), sep="\n")
    
    drugs <- best_labels[i] == drug_tiers
    for (j in 1:nrow(drugs)){
      if (drugs[j,1]){
        text_return <- paste(text_return, sprintf("Drug Tier: %d", as.numeric(drug_tiers[j,2])), sep="\n")
      }
    }
  }
  text_return
}

#########################################
## SHINY UI ############################
#######################################

ui <- dashboardPage(
  dashboardHeader(title = "AMR Dashboard"),
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
                         "Susceptibility",
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
                         "Susceptibility",
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
                               title = "Current Staffing Needs:",
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
                box(width = 12,
                   h3("The Problem and Our Approach"),
                   p("        	Antimicrobial resistance (AMR) occurs when bacteria adapt to administered antibiotics, developing mutations that neutralize or evade their effects. Antimicrobial susceptibility testing is the best option for effective treatment and preventing further resistance; our goal is to provide veterinarians with a tool that leverages data from susceptibility tests to make data-driven decisions for antibody treatment. Therefore, our tool endeavours to be an alternative diagnostic tool that can be referenced before having to use expensive antimicrobial susceptibility testing."),
                   p("To accomplish this, we sought to answer 4 questions: 1) what are the best antibiotic options to use on specific bacteria given the animal species? 2) Does a relationship exist between the age of a studied companion animal and the susceptibility to all the antibiotics tested in the dataset? 3) Are there certain bacteria that are more resistant to antibiotics than others? 4) In the last four years, which antibiotics have domestic canines and felines in New York been most susceptible to? To communicate the findings of our analysis of previous susceptibility tests we endeavored to design a simplistic dashboard that displayed practical information to vets."),
                   h3("Literature Review"),
                   p("Our team first sought to review the literature on antimicrobial resistance (AMR) in general and the One Health approach to AMR so that we could have the background knowledge to approach this project. European surveys have shown that the general public still misunderstands the function and correct use of antibiotics – WHO has highlighted the importance of involving the general public alongside healthcare professionals in combating the emergence of AMR (Zowawi et al., 2015). Antimicrobials are frequently used in food and animal production, and thus AMR has far-reaching effects throughout the entire ecosystem (Kasimanickam et al., 2021). Antibiotic resistance occurs naturally when bacteria are not killed or their growth is not stopped by drugs that were previously effective, and the misuse of antibiotics in humans and animals can accelerate the process (Kasimanickam et al., 2021). Our dashboard seeks to provide vets with the information they need to make data-driven decisions about the appropriate antibiotic to use in a given situation.
"),
                   p("One of our team’s design goals was to make the dashboard as simple and concise as possible so that vets could quickly digest the required information to make a data-driven decision. Naïve Bayesian networks are one of the most effective and simplest Bayesian networks for prediction (Langarizadeh & Maghbeli, 2016). Langarizadeh and Maghbeli (2016) conducted a review of Naïve Bayesian networks and their effectiveness at disease prediction – they found that it had the best performance in most diseases in comparison with other algorithms, and typically had a higher accuracy.
"),
                   p("For our second research question which examines the influence of age on AMR, our team opted to use multiple linear regression instead of binomial logistic regression which has been the favoured method by closely related research ventures. A study by Oguttu et al. (2021) used univariate analysis to assess simple associations between year, season, breed group, age group, sex, and specimen as covariates and extensive drug resistance as the outcome. They found that only the year variable was significantly associated with extensive drug resistance and so it was included in a multivariable logistic model to investigate predictors of extensive drug resistance (Oguttu et al., 2021). This suggests that there may be a temporal factor influencing AMR. Another study sought to model AMR in livestock and companion animals using elastic net logistic regression (Chung et al., 2023). Overall, the study by Chung et al. (2023) also provided us with insight and background information into AMR differences between species, which we did not choose to pursue.
"),
                   h4("References:"),
                   p("Chung, H. C., Foxx, C. L., Hicks, J. A., Stuber, T. P., Friedberg, I., Dorman, K. S., & Harris, B. (2023). An accurate and interpretable model for antimicrobial resistance in pathogenic Escherichia coli from livestock and companion animal species. PLOS ONE, 18(8). https://doi.org/10.1371/journal.pone.0290473
Kasimanickam, V., Kasimanickam, M., & Kasimanickam, R. (2021). Antibiotics use in food animal production: Escalation of antimicrobial resistance: Where are we now in combating amr? Medical Sciences, 9(1), 14. https://doi.org/10.3390/medsci9010014
Langarizadeh, M., & Moghbeli, F. (2016). Applying naive Bayesian networks to disease prediction: A systematic review. Acta Informatica Medica, 24(5), 364. https://doi.org/10.5455/aim.2016.24.364-369
Oguttu, J., Qekwana, D., & Odoi, A. (2021). Prevalence and predictors of antimicrobial resistance among Enterococcus spp.. from dogs presented at a veterinary teaching hospital, South Africa. Frontiers in Veterinary Science, 7. https://doi.org/10.3389/fvets.2020.589439
Zowawi, H. M., Abedalthagafi, M., Mar, F. A., Almalki, T., Kutbi, A. H., Harris-Brown, T., Harbarth, S., Balkhy, H. H., Paterson, D. L., & Hasanain, R. A. (2015). The potential role of social media platforms in community awareness of antibiotic use in the Gulf Cooperation Council states: Luxury or necessity? Journal of Medical Internet Research, 17(10). https://doi.org/10.2196/jmir.3891
"),
                   h3("Description of Methods"),
                   h4("Naïve Bayesian classifier"),
                   p("Our team’s main goal was to help vets make a better decision on what antibiotic to prescribe to pets to reduce AMR spread by giving antibiotics that will be resisted less often. This would mean a lower chance of bacteria becoming resistant to antibiotics rather than when multiple are tested until the correct or susceptible one is found. Thus, the main question our team sought to answer is: what are the best antibiotic options to use on specific bacteria given the animal species? By gathering susceptibility information for each antibiotic, we were able to apply the Naïve Bayesian method; the Naïve Bayesian classifier classifies data based on the probability of it falling into a particular class. Given the species, bacteria, and antibiotic, it calculates the probability that the bacteria is susceptible to the antibiotic and allows the vet to make an educated decision given the data while also considering the drug tiers.
"),
                   h4("Multiple Linear Regression"),
                   p("Another question of importance we strived to answer was: does a relationship exist between the age of a studied companion animal and the susceptibility to all antibiotics tested in the dataset? To answer this question we employed multiple linear regression which leverages explanatory variables to determine the outcome of a dependent variable - our team used the R1…Y1 antibiotic columns in the AMR dataset as the explanatory variables to predict the outcome of age."),
                   h4("Bar Plots"),
                   p("Finally, our team wanted to provide vets with a simple visualization that would allow them to quickly reference the overall resistance and susceptibility rates of each antibiotic and bacteria respectively. To accomplish this, our team created two bar plots – one shows the susceptibility percentage of each antibiotic across all bacteria while the other shows the susceptibility percentage of each bacteria across all antibiotics.
We also attempted to find a trend for which months of the year had the most recorded susceptibility tests for a given bacteria – useful in informing vets when certain bacteria may be most active.
"),
                   )
              )
      )
    )
  )
)

#########################################
## SERVER LOGIC ########################
#######################################
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
      updateActionButton(session, "swapplot_anti", label = "Susceptibility")
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
      plot + ggtitle("Susceptibility Percentage by Antibiotic") +
        xlab("Antibiotic Type") + ylab("Bacteria Susceptibility %") + ylim(0, 100)
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
      updateActionButton(session, "swapplot_bac", label = "Susceptibility")
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
      plot + ggtitle("Susceptibility Percentage by Bacteria") +
        xlab("Bacteria Type") + ylab("Bacteria Susceptibility %") + ylim(0, 100)
    }
  })
  output$best_antibiotics <- renderText({
    if (input$search != "") best_antibiotics(input$search, input$species)
  })
  
  ## LINEAR REGRESSION LOGIC BELOW:
  
  # only grab date of the asked for age
  filter_data <- reactive({
    filtered_data <- tidy_data[
      input$given_age >= tidy_data$estimate - tidy_data$std.error &
        input$given_age <= tidy_data$estimate + tidy_data$std.error, ]
    return(filtered_data)
  })
  
  ###############################
  # DISPLAY AGE ################
  #############################
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
  # DISPLAY TESTS ##############
  #############################
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
}

# Run the application
shinyApp(ui, server)
  
