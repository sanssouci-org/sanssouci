# Define UI for application that draws a histogram
shinyUI(fluidPage(
    tags$head(
        tags$style(HTML('#buttonValidate{background-color:lightblue}'))),  ## this should go to css file?
    
    # Application title
    titlePanel("Permutation-based post hoc confidence bounds for differential gene expression"),
    
    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            splitLayout(
                checkboxInput("checkboxDemo", label = "Use demo data", value = TRUE),
                actionButton("buttonValidate", "Run!"), align = "right"),
            conditionalPanel(condition = "!input.checkboxDemo", 
                             splitLayout(
                                 fileInput("fileData",label = "Input expression data"), 
                                 fileInput("fileCateg", label = "Input condition vector")
                             )
            ), 
            sliderInput("sliderConfLevel", label = "Confidence level", min = 0, 
                            max = 100, value = 90, post = " %"),
            splitLayout(
                selectInput("alternative", label="Alternative", 
                            choices = list("Two sided" = "two.sided", 
                                           "Less"="less", 
                                           "Greater"="greater"),
                            selected = "two.sided"),
                numericInput("numB", label = "Number of permutations", value = 100)
            ),
            splitLayout(
                selectInput("refFamily", label = "Reference family", 
                            choices = list("Simes" = "Simes", "Beta" = "Beta"), 
                            selected = "Simes"),
                uiOutput("inputK")),
            tableOutput("tableBounds"),
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            splitLayout(
                cellWidths = c("60%", "40%"),
                radioButtons("choiceYaxis", label = "'y' axis label", 
                             choices = list("p-values" = "pval", 
                                            "Adjusted p-values" = "adjPval",
                                            "Calibration thresholds" = "thr"), 
                             selected = "pval",
                             inline = TRUE),
                checkboxInput("symetric", 
                              label = "Symetric foldchange threshold", 
                              value = TRUE)),
            plotly::plotlyOutput("volcanoplot", height = "600px"), 
            fluidRow(
                actionButton("resetCSV", "Reset Selections"), 
                downloadButton("downloadData", "Download binary csv of user selection")
            )
        )
    ),
    p(em("This interactive ",
      a("shiny", href = "https://shiny.rstudio.com"),
      "application is developed by",
      "Nicolas Enjalbert-Courrech", 
      "and",
      a("Pierre Neuvial", href = "https://www.math.univ-toulouse.fr/~pneuvial/"),
      "for the R package ",
      a("sansSouci.", href = "https://pneuvial.github.io/sanssouci/"),
      "It implements permutation-based post hoc inference bounds for differential gene expression analysis, see dedicated ",
      a("vignette.", href = "https://pneuvial.github.io/sanssouci/articles/post-hoc_differential-expression.html"), 
      "The source code for this app is currently available from ",
      a("this url.", href = "https://github.com/pneuvial/sanssouci/tree/volcano-plot-app/inst/shiny-examples/volcano-plot"), 
      "For any question, please file an",
      a("issue.", href = "https://github.com/pneuvial/sanssouci/issues"))
)))
