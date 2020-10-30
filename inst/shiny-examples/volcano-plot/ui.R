library(shiny)
library(plotly)
library(sansSouci)
library(sansSouci.data)
library(htmlwidgets)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("Post hoc confidence bounds for volcano plots"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            checkboxInput("checkboxDemo", label = "Use demo data", value = FALSE),
            conditionalPanel(condition = "!input.checkboxDemo", 
                             fileInput("fileData",label = "Input data"), 
                             fileInput("fileCateg", label = "Input categ")
            ), 
            sliderInput("sliderAlpha", label = "alpha", min = 0, 
                        max = 1, value = 0.05, step = 0.01),
            
            # actionButton("buttonValidateAlpha", "Go - Validation of alpha"),
            tags$p("Chosen threshold values :"), 
            textOutput("thresholdTxt")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            checkboxInput("symetric", label = "Symetric foldchange thresold", value = TRUE),
            tableOutput("tableBounds"),
            # tableOutput("datatable"), 
            plotlyOutput("volcanoplot")
            # textOutput("text1"), 
            # textOutput("text2"), 
            # textOutput("text12")
        )
    )
))
