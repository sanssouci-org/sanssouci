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
# data(hgu95av2_GO_BP, package = "sansSouci.data")
# data(hgu95av2_GO_MF, package = "sansSouci.data")
# data(hgu95av2_GO_CC, package = "sansSouci.data")
data(expr_ALL_GO, package = "sansSouci.data")

shinyUI(fluidPage(
  
  includeCSS("www/style.css"),
  
  # Application title
  titlePanel("IIDEA: Interactive Inference for Differential Expression Analyses"),
  # Sidebar with panel
  sidebarLayout(
    sidebarPanel(
      # wellPanel(
      splitLayout(
        htmlOutput("help"),
        checkboxInput("checkboxDemo", 
                      label = "Use demo data", 
                      value = TRUE),
        # uiOutput("CheckData"),
        actionButton("buttonValidate", "Run!", )),
      conditionalPanel(condition = "!input.checkboxDemo", 
                       fileInput("fileData",
                                 label = p("Gene expression data matrix", 
                                           bsButton("QfileData", 
                                                    label = "", 
                                                    icon = icon("question"), 
                                                    style = "info", 
                                                    size = "extra-small"))),
                       bsTooltip("QfileData", "Upload a CSV file containing matrix with genes in rows and samples in column. Column names should be in (in {0, 1})",
                                 "right", 
                                 options = list(container = "body"), 
                                 trigger = "hover")
      ), 
      conditionalPanel(condition = "!input.checkboxDemo",
                       # splitLayout(
                       #   fileInput("fileAnnotation", label = p("Input gene annotation",  
                       #                                         bsButton("QfileAnnotation", label = "", icon = icon("question"), style = "info", size = "extra-small"))), 
                       fileInput("fileGroup", 
                                 label = p("Gene set matrix",
                                           bsButton("QfileGroup", 
                                                    label = "", 
                                                    icon = icon("question"), 
                                                    style = "info", 
                                                    size = "extra-small"))),
                       # ),
                       # bsTooltip(id = "QfileAnnotation", title = 'Upload a CSV file containing matrix within two columns. One called "Id" contains index label from matrix and the other, called "nameGene", contains names of associated genes.', placement = "bottom",  options = list(container = "body"), trigger = "focus"),
                       bsTooltip(id = "QfileGroup", 
                                 title = 'Upload a RDS file containing matrix within nameGenes in line index. Binary vector composed this matrix for each gene set.', 
                                 placement = "bottom", 
                                 trigger = "hover", 
                                 options = NULL),
                       downloadButton("downloadExampleData", "Download example data set")
      ), 
      sliderInput("sliderConfLevel", 
                  "Confidence level", 
                  # label = p("Confidence level", 
                  #           bsButton("QconfLevel", 
                  #                    label = "", 
                  #                    icon = icon("question"), 
                  #                    style = "info", 
                  #                    size = "extra-small")),
                  min = 0, 
                  max = 100, value = 90, post = " %"),
      # bsTooltip("QconfLevel", "Confidence level",
      #           "right", 
      #           options = list(container = "body"), 
      #           trigger = "hover"),
      checkboxInput("checkboxAdvancedParam", 
                    label = p("Advanced parameters"),
                              # bsButton("Qparam", 
                              #          label = "", 
                              #          icon = icon("question"), 
                              #          style = "info", 
                              #          size = "extra-small")),
                    value = FALSE),
      # bsTooltip(id = "Qparam", 
      #           title = paste("Select parameters to implement permutation-based post hoc inference bounds for differential gene expression analysis, see dedicated ", 
      #                           a("vignette.", 
      #                             href = "https://pneuvial.github.io/sanssouci/articles/post-hoc_differential-expression.html")), 
      #           trigger = c("click", "hover"),
      #           options = NULL),
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
                   
      ),
      # )
    ),
    
    # Main panel
    mainPanel(
      flowLayout(
        selectInput("choiceYaxis", label = "'y' axis label", 
                  choices = list("p-values" = "pval", 
                                 "Adjusted p-values" = "adjPval",
                                 "Number of false positves" = "thr"), 
                  selected = "pval"),
      checkboxInput("symetric", 
                    label = "Symmetric fold change threshold", 
                    value = FALSE),
      h2("Volcano plot",
         bsButton("Qparam1", label = "", icon = icon("question"), style = "info", size = "extra-small")),
      bsPopover(id = "Qparam1", 
                title = "VolcanoPlot", 
                content = paste('Select genes by dragging horizontal or vertical bars, of using "box select" or "lasso select" from the plot menu. The table in the left panel gives post-hoc bounds for these selections.'), 
                placement = "bottom", 
                trigger = "hover", 
                options = NULL)),
      conditionalPanel(condition = "input.tabSelected==1",
                         plotly::plotlyOutput("volcanoplotPosteriori", height = "600px"), 
                       
                       fluidRow(
                         actionButton("resetCSV", "Reset Selections"), 
                         downloadButton("downloadData", "Download csv file with user selection")
                         # bsButton("Qdownload", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                         # bsTooltip("Qdownload", "Delete your select manual selection from the post hoc bounds and downloadable csv file. Download a csv file containing matrix with binary vector of your User selection",
                         #           "right", options = list(container = "body"), trigger = "focus"),
                         # bsTooltip(id = "resetCSV", title = "Delete you box select manual selection from the post hoc bounds and downloadable csv file.", placement = "bottom", trigger = "hover", options = NULL),
                         # bsTooltip(id = "downloadData", title = "Download a csv file containing matrix with binary vector of your User selection", placement = "bottom", trigger = "hover", options = NULL)
                      )),
                       
                       # splitLayout(
                       #   plotlyOutput("curveMaxFPBoth"), 
                       #   plotlyOutput("curveMaxFPSelect")
                       # )), 
                       # 
      conditionalPanel(condition = "input.tabSelected==2",
                       plotly::plotlyOutput("volcanoplotPriori", height = "600px")
                       # plotlyOutput(outputId = "curveMaxFPGroup")
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
