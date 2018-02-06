library("plotly")
library("sansSouci")  ## devtools::install_github("pneuvial/sanssouci@develop")

data(volcano, package = "sansSouci")
dataSets <- unique(volcano[["dataSet"]])
volcano[["logp"]] <- -log10(volcano[["p.value"]])

ui <- fluidPage(
    inputPanel(
        selectInput("dataSet", "Data set", choices = dataSets, selected = "golub"),
        numericInput("alpha", "Target confidence level:", 0.1, min = 0, max = 1)),
    wellPanel(h3(textOutput("bound")),
              plotlyOutput("plot"))
)

server <- function(input, output, session) {
    
    volc <- reactive({ 
        volcano[which(volcano$dataSet == input$dataSet), ]
    })
    
    output$plot <- renderPlotly({
        datly <- volc()
        rg <- range(datly$meanDiff)
        rg <- max(abs(rg))*c(-1,1)
        print(rg)
        
        #Vbar0 <- upperBoundFP(sort(x, decreasing=TRUE), thr)
        plot_ly(datly, type = "scatter", mode = "markers",
                x = ~meanDiff, y = ~logp, key = ~id, 
                text = ~id, hoverinfo = 'text') %>%
            layout(dragmode = "select",
                   xaxis = list(range = rg))
    })
    
    output$bound <- renderText({
        datly <- volc()
        alpha <- input$alpha
        d <- event_data("plotly_selected")
        msg <- "Select a set of points"
        if (is.null(d)) {
            mm <- which(datly[["adjp"]] <= alpha)
        } else {
            mm <- match(d$key, datly[["id"]])
            Vbar <- posthocBySimes(datly[["p.value"]], mm, alpha)
            msg <- sprintf("At least %s true positives among %s selected genes (FDP <  %s)", Vbar, nrow(d), round(1-Vbar/nrow(d), 2))
        }
        #msg <- paste(msg, str(d), collapse = "\n")
        msg
    })
}
## TODO: 
## x select data set
## * export selected genes
## * multiple selection: exists since dec 17 in js, not in R yet
shinyApp(ui, server)