ui <- fluidPage(
    titlePanel("Post hoc confidence bounds for volcano plots"),
    inputPanel(
        selectInput("dataSet", "Data set", choices = dataSets, selected = "chiaretti"),
        numericInput("alpha", "Target confidence level:", 0.05, min = 0, max = 1)),
    fluidRow(
        column(9,
               wellPanel(h3(textOutput("bound")),
                         plotlyOutput("plot"))),
        column(3,
               wellPanel(h3("selected genes"),
                         textOutput("brush")),
               wellPanel(h3("clicked genes"),
                         textOutput("click"))
        )
    )
)

server <- function(input, output, session) {
    
    volc <- reactive({ 
        dat <- volcano[which(volcano$dataSet == input$dataSet), ]
        subset(dat, abs(meanDiff) > 0.2 | logp < -1)
    })
    
    output$plot <- renderPlotly({
        datly <- volc()
        rg <- range(datly$meanDiff)
        rg <- max(abs(rg))*c(-1,1)

        d <- event_data("plotly_selected")
        print(d)
        mm <- match(d$key, datly[["id"]])
        print(mm)
        datly$selected <- 0
        datly$selected[mm] <- 1
        table(datly$selected)
        
        library("ggplot2")
        p <- ggplot(datly, aes(x = meanDiff, y = -logp, 
                               colour = factor(selected), key = id))
        p <- p + geom_point(alpha = 0.2)
        p <- p + xlim(range(datly$meanDiff))
        p <- p + ylim(c(0, max(-datly$logp)))
        ggplotly(p)  %>% 
            config(displayModeBar = F) %>%
            layout(dragmode = "select",
                   xaxis = list(range = rg))
    })

    output$brush <- renderText({
        d <- event_data("plotly_selected")
        if (is.null(d)) {
            "Nothing to show here yet. Need a selection."
        } else {
            d$key
        }
    })
    
    output$bound <- renderText({
        datly <- volc()
        alpha <- input$alpha
        d <- event_data("plotly_selected")
        msg <- "Select a set of points"
        if (is.null(d)) {
            ## TODO: default selection
            ## mm <- which(datly[["adjp"]] <= alpha)
        } else {
            mm <- match(d$key, datly[["id"]])
            if (!all(is.na(mm))) {
                Vbar <- posthocBySimes(datly[["p.value"]], mm, alpha)
                msg <- sprintf("At least %s true positives among %s selected genes", Vbar, nrow(d))
                fdp <- round(1 - Vbar/nrow(d), 2)
                msg <- sprintf("%s (FDP %s %s)", msg, ifelse(fdp==0, "=", "<="), fdp)
            }
        }
        # msg <- paste(msg, str(d), collapse = "\n")
        msg
    })
    
    output$click <- renderPrint({
        d <- event_data("plotly_click")
        if (is.null(d)) "Click events appear here (double-click to clear)" else d
    })
}
## TODO: 
## x select data set
## * export selected genes
## * multiple selection: exists since dec 17 in js, not in R yet
## * on click: afficher les termes du GO correspondants
#shinyApp(ui, server)