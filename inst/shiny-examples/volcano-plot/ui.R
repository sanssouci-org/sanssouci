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
                             splitLayout(
                                 fileInput("fileData",label = "Input data"), 
                                 fileInput("fileCateg", label = "Input categ")
                             )
            ), 
            splitLayout(
                sliderInput("sliderAlpha", label = "alpha", min = 0, 
                            max = 1, value = 0.05, step = 0.01),
                numericInput("numB", label = "Number of permutation", value = 100)
            ), 
            splitLayout(
                selectInput("refFamily", label = "RefFamily", 
                            choices = list("Simes" = "Simes", "Beta" = "Beta"), 
                            selected = "Simes"),
                selectInput("alternative", label="Alternative", 
                            choices = list("Two sided"= "two.sided", "Less"="less", "Greater"="greater"),
                            selected="two.sided")
            ),
            uiOutput("inputK"),
            
            actionButton("buttonValidate", "Go - Validation of inputs"),
            tableOutput("tableBounds")
            
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            splitLayout(
                checkboxInput("symetric", label = "Symetric foldchange thresold", value = TRUE),
                radioButtons("choiceYaxis", label = "Y axis", 
                             choices = list("P values"="pval", 
                                            "Adjusted p values"="adjPval",
                                            "Calibration threshold"="thr"), 
                             selected="pval"
                )
            ),
            plotlyOutput("volcanoplot", height = "600px"), 
            fluidRow(
                actionButton("resetCSV", "Reset Selections"), 
                downloadButton("downloadData", "Download binary csv of user selection")
            )
        )
    )
))
