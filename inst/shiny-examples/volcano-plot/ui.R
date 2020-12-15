library(shiny)
library(plotly)
library(sansSouci)
library(sansSouci.data)  
stopifnot(packageVersion("sansSouci.data") >= '0.2.2')
library(ggplot2)
library(dplyr)
library(htmlwidgets)
library(DT)

data(expr_ALL, package = "sansSouci.data")
data(expr_ALL_annotation, package = "sansSouci.data")
data(hgu95av2_GO_BP, package = "sansSouci.data")
# data(hgu95av2_GO_MF, package = "sansSouci.data")
# data(hgu95av2_GO_CC, package = "sansSouci.data")

shinyUI(fluidPage(
  
  includeCSS("www/style.css"),
  
  # Application title
  titlePanel("Permutation-based post hoc confidence bounds for differential gene expression"),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      # wellPanel(
      splitLayout(
        checkboxInput("checkboxDemo", label = "Use demo data", value = TRUE),
        actionButton("buttonValidate", "Perform calibration!"), align = "right"),
      conditionalPanel(condition = "!input.checkboxDemo", 
                       splitLayout(
                         fileInput("fileData",label = "Input expression data"), 
                         fileInput("fileCateg", label = "Input condition vector")
                       )
      ), 
      conditionalPanel(condition = "!input.checkboxDemo",
                       splitLayout(
                         fileInput("fileAnnotation", label = "Input gene annotation"), 
                         fileInput("fileGroup", label = "Matrix of biological function")
                       )
      ), 
      sliderInput("sliderConfLevel", label = "Confidence level", min = 0, 
                  max = 100, value = 90, post = " %"),
      splitLayout(
        selectInput("alternative", label = "Alternative", 
                    choices = list("Two sided" = "two.sided", 
                                   "Less" = "less", 
                                   "Greater" = "greater"),
                    selected = "two.sided"),
        numericInput("numB", label = "Number of permutations", value = 100, min = 2)
      ),
      splitLayout(
        selectInput("refFamily", label = "Reference family", 
                    choices = list("Simes" = "Simes", "Beta" = "Beta"), 
                    selected = "Simes"),
        uiOutput("inputK")#)
      ),
      tabsetPanel( id = "tabSelected",
        tabPanel("User selections", value = 1,
                 DTOutput("tableBounds")),
        tabPanel("Gene sets", value = 2,
                 DTOutput("tableBoundsGroup") )
        
      )
      # )
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
                      label = "Symmetric fold change threshold", 
                      value = FALSE)), 
      conditionalPanel(condition = "input.tabSelected==2",
                       uiOutput("choiceGroupUI")),
      
      conditionalPanel(condition = "input.tabSelected==1",
                       plotly::plotlyOutput("volcanoplotPosteriori", height = "600px"), 
                       fluidRow(
                         actionButton("resetCSV", "Reset Selections"), 
                         downloadButton("downloadData", "Download binary csv of user selection")
                       ), 
                       splitLayout(
                         plotlyOutput("curveMaxFPBoth"), 
                         plotlyOutput("curveMaxFPSelect")
                       )), 
      conditionalPanel(condition = "input.tabSelected==2",
                       plotly::plotlyOutput("volcanoplotPriori", height = "600px"), 
                       plotlyOutput(outputId = "curveMaxFPGroup")
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
       "The source code for this app is available from ",
       a("this url.", href = "https://github.com/pneuvial/sanssouci/tree/develop/inst/shiny-examples/volcano-plot"), 
       "For any question, please file an",
       a("issue.", href = "https://github.com/pneuvial/sanssouci/issues"))
  )))
