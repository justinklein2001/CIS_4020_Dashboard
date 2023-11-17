#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

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

drug_tiers <- read_csv(drug_tiers_data_file)

parquet_file <- arrow::read_parquet(amr_data_file)
print(unique(parquet_file$R1))

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

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  title = "Dashboard",
  
  # Sidebar with a slider input for number of bins 
  fluidRow(
    column(width = 6, offset = 0, align = "center",
           style = "background-color:#d3d3d3;border: 3px solid #000000;",
           plotOutput("antiPlot")
    ),
    column(width = 6, offset = 0, align = "center",
           style = "background-color:#d3d3d3;border: 3px solid #000000;",
           plotOutput("bacPlot")
    ),
  ),  
  fluidRow(
    column(width = 4, offset = 0, align = "center", 
           style = "background-color:#d3d3d3;border: 3px solid #000000; height: 10%",
           #plotOutput("antiPlot"),
           sliderInput("bins",
                       "# of Antibiotics:",
                       min = 4,
                       max = nrow(res_prob),
                       value = 20,
                       step = 1,
                       width = "75%"),
    ),
    column(width = 2, offset = 0,
           style = "background-color:#d3d3d3;border: 3px solid #000000;height: 105.5px;",
           align = "center",
           actionButton(
             "swapplot_anti", 
             "Subsceptibility",
             style='padding:4px; font-size:80%; margin-top: 15%;'
           ),
           actionButton(
             "sorting_anti", 
             "Ascending",
             style='padding:4px; font-size:80%; margin-top: 15%;'
           )
    ),
    column(width = 4, offset = 0, align = "center", 
           style = "background-color:#d3d3d3;border: 3px solid #000000; height: 10%",
           #plotOutput("antiPlot"),
           sliderInput("bins_bac",
                       "# of Bacteria:",
                       min = 2,
                       max = 6,#nrow(res_collection) - 1,
                       value = 4,
                       step = 1,
                       width = "75%"),
    ),
    column(width = 2, offset = 0,
           style = "background-color:#d3d3d3;border: 3px solid #000000;height: 105.5px;",
           align = "center",
           actionButton(
             "swapplot_bac", 
             "Subsceptibility",
             style='padding:4px; font-size:80%; margin-top: 15%;'
           ),
           actionButton(
             "sorting_bac", 
             "Ascending",
             style='padding:4px; font-size:80%; margin-top: 15%;'
           )
    ),
  
  # Show a plot of the generated distribution
  #mainPanel(
    #plotOutput("antiPlot")
  #)
  ),
  fluidRow(align = "center", style = "background-color:#d3d3d3;",
    column(width = 3, offset = 4, align = "center",
      selectizeInput(
        inputId = "search", 
        label = "Best Antibiotics Search",
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
    column(width = 1, offset = 0, align = "center",
      selectInput(
        inputId = "species",
        label = "Species",
        choices = c("Canine", "Feline"),
        selected = "Canine",
        multiple = FALSE,
        selectize = TRUE,
      ),
      textOutput("text"),
    )
  ),
  fluidRow(align = "center", style = "background-color:#d3d3d3;",
    verbatimTextOutput("best_antibiotics"),
  ),
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
                                           y=chance)) + geom_bar(stat = "identity")
    }
    else{
      plot <- ggplot(bins, aes(x=reorder(rownames(bins), chance), 
                               y=chance)) + geom_bar(stat = "identity")
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
                               y=prob)) + geom_bar(stat = "identity")
    }
    else{
      plot <- ggplot(bins_bac, aes(x=reorder(org_standard, prob), 
                               y=prob)) + geom_bar(stat = "identity")
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

}

# Run the application 
shinyApp(ui = ui, server = server)
