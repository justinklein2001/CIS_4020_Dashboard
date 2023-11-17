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
colnames(res_prob) <- 'chance'

sub_prob <- as.data.frame(sub_prob)
colnames(sub_prob) <- 'chance'

#res_prob <- arrange(res_prob, desc(resisted))
#print(res_prob)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  title = "Dashboard",
  
  # Sidebar with a slider input for number of bins 
  fluidRow(
    column(width = 4, offset = 0, align = "left",
           style = "background-color:#d3d3d3;border: 3px solid #000000;",
           plotOutput("antiPlot")
    ),
  ),  
  fluidRow(
    column(width = 3, offset = 0, align = "center", 
           style = "background-color:#d3d3d3;border: 3px solid #000000;",
           #plotOutput("antiPlot"),
           sliderInput("bins",
                       "# of Antibiotics:",
                       min = 5,
                       max = nrow(res_prob) - 1,
                       value = 20,
                       step = 1,
                       width = "75%"),
    ),
    column(width = 1, offset = 0,
           style = "background-color:#d3d3d3;border: 3px solid #000000;height: 105.5px;",
           align = "center",
           actionButton(
             "swapplot", 
             "Subsceptibility",
             style='padding:4px; font-size:80%; margin-top: 27%;'
           ),
           actionButton(
             "sorting", 
             "Ascending",
             style='padding:4px; font-size:80%;'
           )
    ),
  
  # Show a plot of the generated distribution
  #mainPanel(
    #plotOutput("antiPlot")
  #)
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  params <- reactiveValues(
    prob = arrange(res_prob, desc(chance)), 
    descending = TRUE, 
    res = TRUE
  )
  
  observeEvent(input$sorting, {
    if (params$descending){
      params$prob <- arrange(params$prob, chance)
      params$descending <- FALSE
      updateActionButton(session, "sorting", label = "Descending")
    }
    else{
      params$prob <- arrange(params$prob, desc(chance))
      params$descending <- TRUE
      updateActionButton(session, "sorting", label = "Ascending")
    }
  })
  
  observeEvent(input$swapplot, {
    updateActionButton(session, "swapplot", label = "Hello")
    if (params$res){
      params$prob <- sub_prob
      params$res <- FALSE
      updateActionButton(session, "swapplot", label = "Subsceptibility")
    }
    else{
      params$prob <- res_prob
      params$res <- TRUE
      updateActionButton(session, "swapplot", label = "Resistance")
    }
  })
  
  output$antiPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    bins <- as.data.frame(params$prob[0:input$bins + 1,, drop = FALSE])
    
    if (params$descending){
      plot <- ggplot(bins, aes(x=reorder(rownames(bins), -chance), 
                                           y=chance)) + geom_bar(stat = "identity")
    }
    else{
      plot <- ggplot(bins, aes(x=reorder(rownames(bins), chance), 
                               y=chance)) + geom_bar(stat = "identity")
    }
    
    if (params$res){
      plot + ggtitle("Resistance Percentage by Antibiotic") +
        xlab("Antibiotic Type") + ylab("Bacteria Resistance %") + ylim(0, 100)
    }
    else{
      plot + ggtitle("Subsceptibility Percentage by Antibiotic") +
        xlab("Antibiotic Type") + ylab("Bacteria Subsceptibility %") + ylim(0, 100)
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)