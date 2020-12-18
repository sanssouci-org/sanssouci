library(shiny)
library(plotly)
library(sansSouci)
library(sansSouci.data)  
stopifnot(packageVersion("sansSouci.data") >= '0.2.2')
library(ggplot2)
library(dplyr)
library(htmlwidgets)
library(DT)
library(shinyBS)
library(stringr)

data(expr_ALL, package = "sansSouci.data")
#data(expr_ALL_annotation, package = "sansSouci.data")
data(hgu95av2_GO_BP, package = "sansSouci.data")
# data(hgu95av2_GO_MF, package = "sansSouci.data")
# data(hgu95av2_GO_CC, package = "sansSouci.data")
data(expr_ALL_GO, package = "sansSouci.data")

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
        # uiOutput("CheckData"),
        actionButton("buttonValidate", "Perform calibration!"),
        align = "right"),
      conditionalPanel(condition = "!input.checkboxDemo", 
                       fileInput("fileData",label = p("Input expression data", 
                                                      bsButton("QfileData", label = "", icon = icon("question"), style = "info", size = "extra-small"))),
                       bsTooltip("QfileData", "Upload a RDS file containing matrix within  __ in line index and categories in column (in {0, 1})",
                                 "right", options = list(container = "body"), trigger = "focus")
      ), 
      conditionalPanel(condition = "!input.checkboxDemo",
                       # splitLayout(
                       #   fileInput("fileAnnotation", label = p("Input gene annotation",  
                       #                                         bsButton("QfileAnnotation", label = "", icon = icon("question"), style = "info", size = "extra-small"))), 
                       fileInput("fileGroup", label = p("Matrix of biological function",
                                                        bsButton("QfileGroup", label = "", icon = icon("question"), style = "info", size = "extra-small"))),
                       # ),
                       # bsTooltip(id = "QfileAnnotation", title = 'Upload a RDS file containing matrix within two columns. One called "Id" contains index label from matrix and the other, called "nameGene", contains names of associated genes.', placement = "bottom",  options = list(container = "body"), trigger = "focus"),
                       bsTooltip(id = "QfileGroup", title = 'Upload a RDS file containing matrix within nameGenes in line index. Binary vecotr composed this matrix for each biological function.', placement = "bottom", trigger = "hover", options = NULL)
      ), 
      bsButton("Qparam", label = "", icon = icon("question"), style = "info", size = "extra-small"),
      bsPopover(id = "Qparam", title = "Parameters", content = paste("Select parameters to implement permutation-based post hoc inference bounds for differential gene expression analysis, see dedicated ", a("vignette.", href = "https://pneuvial.github.io/sanssouci/articles/post-hoc_differential-expression.html")), 
                placement = "bottom", trigger = "focus", options = NULL),
      
      checkboxInput("checkboxAdvancedParam", label = "Advanced calibration parameters", value = FALSE),
      sliderInput("sliderConfLevel", label = "Confidence level", min = 0, 
                  max = 100, value = 90, post = " %"),
      conditionalPanel(condition = "input.checkboxAdvancedParam",
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
                       )),
      tabsetPanel( id = "tabSelected",
                   tabPanel("User selections", value = 1,
                            
                            uiOutput("OutQtableBounds"),
                            DTOutput("tableBounds")),
                   tabPanel("Gene sets", value = 2,
                            uiOutput("OutQtableBoundsGroup"),
                            DTOutput("tableBoundsGroup") )
                   
      )
      # )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      splitLayout(
        cellWidths = c("55%", "35%", "10%"),
        radioButtons("choiceYaxis", label = "'y' axis label", 
                     choices = list("p-values" = "pval", 
                                    "Adjusted p-values" = "adjPval",
                                    "Calibration thresholds" = "thr"), 
                     selected = "pval",
                     inline = TRUE),
        
        checkboxInput("symetric", 
                      label = "Symmetric fold change threshold", 
                      value = FALSE),
        bsButton("Qparam1", label = "", icon = icon("question"), style = "info", size = "extra-small"),
        bsPopover(id = "Qparam1", title = "VolcanoPlot", content = paste('You can make some selections on the volcano plot. Firstly, you can draggle threshold (orange pointed lines) to select both upper right and left corner. Post-hoc bounds of these selections are printed in the table on the left of the app. Secondly, you can select some point with the "box select" or the "lasso select". Post-hoc bounds are also on the left of the app.'), 
                  placement = "bottom", trigger = "focus", options = NULL)
        
      ),
      
      
      conditionalPanel(condition = "input.tabSelected==1",
                       plotly::plotlyOutput("volcanoplotPosteriori", height = "600px"), 
                       fluidRow(
                         actionButton("resetCSV", "Reset Selections"), 
                         downloadButton("downloadData", "Download binary csv of user selection"), 
                         bsButton("Qdownload", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                         bsTooltip("Qdownload", "Delete you box select manual selection from the post hoc bounds and downloadable csv file. Download a csv file containing matrix with binary vector of your User selection",
                                   "right", options = list(container = "body"), trigger = "focus"),
                         bsTooltip(id = "resetCSV", title = "Delete you box select manual selection from the post hoc bounds and downloadable csv file.", placement = "bottom", trigger = "hover", options = NULL),
                         bsTooltip(id = "downloadData", title = "Download a csv file containing matrix with binary vector of your User selection", placement = "bottom", trigger = "hover", options = NULL)
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
