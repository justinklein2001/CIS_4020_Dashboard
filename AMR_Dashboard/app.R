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


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 50,
                        value = 30)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("monthTestPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$monthTestPlot <- renderPlot({
      # subject to change based off final implementation
      amr_data_file <- "./data/AMR Data.parquet"
      
      
      # Read the Parquet file
      parquet_file <- arrow::read_parquet(amr_data_file)
      
      # store counts of each type of year in a table
      test_counts <- table(parquet_file[["order_month"]])
      
      months_count_df <- data.frame(value = as.numeric(names(test_counts)), count = as.numeric(test_counts))
      
      # plot that baby
      ggplot(months_count_df, aes(x = value, y = count)) +
        geom_line() +
        geom_point() +
        labs("Number of Tests Per Month",
             x = "Months",
             y = "# of Tests")
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
