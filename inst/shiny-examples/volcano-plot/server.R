library(shiny)
library(plotly)
library(sansSouci)
library(sansSouci.data)
data(expr_ALL, package = "sansSouci.data")
data(expr_ALL_annotation, package = "sansSouci.data")

shinyServer(function(input, output) {
    
    library(sansSouci)
    library(ggplot2)
    library(plotly)
    
    source("function.R")
    
    options(shiny.maxRequestSize=1024^3)
    
    data <- reactive (
        # eventReactive(input$buttonValidateAlpha, 
        {
            data <- list()
            if(input$checkboxDemo){
                # data_url <- url("https://plmbox.math.cnrs.fr/f/755496cc4c154a6dbab0/?dl=1")
                # data$matrix <- read(data_url)
                # data(expr_ALL, package = "sansSouci.data")
                data$matrix <- expr_ALL
                
                
                categ <- colnames(data$matrix)
                data$categ <- rep(0, length(categ))
                data$categ[which(categ == "NEG")] <- 1
                
                data$annotation <- expr_ALL_annotation[c('affy_hg_u95av2','hgnc_symbol')]
                
            } else {
                req(input$fileData)
                file <- req(input$fileData)
                data$matrix <- readRDS(file$datapath)
                
                req(input$fileCateg)
                fileCateg <- req(input$fileCateg)
                data$categ <- readRDS(fileCateg$datapath)
            }
            return(data)
        }
    )
    output$test <- renderPrint({input$checkboxDemo})
    output$datatable <- renderTable({head(data()$matrix)})
    
    alpha <- reactive(
        # eventReactive(
        # input$buttonValidateAlpha,
        {
            req(input$sliderAlpha)
        })
    cal <- reactive({
        calibrateJER(data()$matrix, categ = data()$categ, B = 1e2, alpha = alpha(), refFamily="Simes")
    })
    
    df <- reactive({
        dex <- rowWelchTests(data()$matrix, data()$categ)
        pval <- dex[["p.value"]]
        logp <- -log10(pval)
        fc <- dex$meanDiff
        adjp <- p.adjust(pval, method = "BH")
        return(list(logp=logp, fc=fc, adjp=adjp, pval=pval))
    })
    
    
    
    
    # volcanoplotReactive <- eventReactive(input$buttonVolcano, {
    #     req(input$sliderAdjPvalue)
    #     req(input$sliderFoldChange)
    #     VolcanoPlotNico(data()$matrix, thr = cal()$thr, categ = data()$categ, q=req(input$sliderAdjPvalue), r=req(input$sliderFoldChange))
    #     
    # })
    
    
    vertical <- reactive({event_data("plotly_relayout", source='A')})
    
    xint <- reactiveVal(0.5)
    observeEvent(vertical()[["shapes[0].x0"]], {
        if(input$symetric){
            newValue <-  vertical()[["shapes[0].x0"]]
            xint(newValue)
            xint2(-newValue)
        } else {
            newValue <-  vertical()[["shapes[0].x0"]]   
            xint(newValue)  
        }
    })
    xint2 <- reactiveVal(-0.5)
    observeEvent(vertical()[["shapes[2].x0"]], {
        if( input$symetric){
            newValue <-  vertical()[["shapes[2].x0"]]
            xint( - newValue)
            xint2(newValue)
        } else {
            newValue <- vertical()[["shapes[2].x0"]] 
            xint2(newValue)
        }
        
    })
    
    yint <- reactiveVal(-log10(0.05))
    observeEvent(vertical()[["shapes[1].y0"]], {
        p <- 10^(-vertical()[["shapes[1].y0"]])
        y_sel <- which((df()$pval <= p))          ## selected by  p-value
        newValue <- Inf
        if (length(y_sel) > 0) {
            newValue <- min(df()$logp[y_sel])              ## threshold on the log(p-value) scale
        }
        yint(newValue)
    })
    
    output$thresholdTxt <- renderText({
        paste(c(sprintf("Foldchange:  right %s left %s . P_values %s", 
                        round(xint(), digits = 3), round(xint2(), digits = 3), round(10^(-yint()), digits = 3)
        )))
    })
    
    selectedGenes <- reactive({
        req(df())
        ## gene selections
        sel1 <- which(df()$logp >= yint() & df()$fc >= xint()) 
        sel2 <- which(df()$logp >= yint() & df()$fc <= xint2())
        sel12 <- sort(union(sel1,sel2))
        return(list(sel1=sel1, sel2=sel2, sel12=sel12))
    })
    
    TP_FDP <- reactive({
        req(selectedGenes())
        ## post hoc bounds in selections
        n1 <- length(selectedGenes()$sel1)
        FP1 <- maxFP(df()$pval[selectedGenes()$sel1], thr = cal()$thr)
        TP1 <- n1 - FP1
        FDP1 <- round(FP1/max(n1, 1), 2)
        
        n2 <- length(selectedGenes()$sel2)
        FP2 <- maxFP(df()$pval[selectedGenes()$sel2], thr = cal()$thr)
        TP2 <- n2 - FP2
        FDP2 <- round(FP2/max(n2, 1), 2)
        
        n12 <- length(selectedGenes()$sel12)
        FP12 <- maxFP(df()$pval[selectedGenes()$sel12], thr = cal()$thr)
        TP12 <- n12 - FP12
        FDP12 <- round(FP12/max(n12, 1), 2)
        
        return(list(n1=n1, TP1=TP1, FDP1=FDP1, n2=n2, TP2=TP2, FDP2=FDP2, n12=n12, TP12=TP12, FDP12=FDP12))
        
    })
    
    selected_points <- reactive({
        list(x = df()$fc[selectedGenes()$sel12], 
             y = df()$logp[selectedGenes()$sel12])
    })
    
    output$tableBounds <- renderTable({
        req(TP_FDP())
        tab <- data.frame(`Selection`=c("Upper Right", "Upper Left", "Both parts"), 
                          "Number of genes"=c(TP_FDP()$n1, TP_FDP()$n2, TP_FDP()$n12),
                          `Lower bound on true positives`=as.integer(c(TP_FDP()$TP1, TP_FDP()$TP2, TP_FDP()$TP12)), 
                          `Upper bound on false discovery proportion` = c(TP_FDP()$FDP1, TP_FDP()$FDP2, TP_FDP()$FDP12))
        return(tab)
    })
    
    # output$text1 <- renderText({
    #     req(TP_FDP())
    #     lte <- "≤"
    #     gte <- "≥"
    #     paste(c(sprintf("Up Right: %s genes\nTP %s %s\nFDP %s %s", req(TP_FDP()$n1), gte,
    #                     req(TP_FDP()$TP1), lte, req(TP_FDP()$FDP1))))
    # })
    # output$text12 <- renderText({
    #     req(TP_FDP())
    #     lte <- "≤"
    #     gte <- "≥"
    #     paste(c(sprintf("Both parts: %s genes selected\nAt least %s true positives (FDP %s %s)",
    #                     req(TP_FDP()$n12), req(TP_FDP()$TP12), lte, req(TP_FDP()$FDP12))))
    # })
    # output$text2 <- renderText({
    #     req(TP_FDP())
    #     lte <- "≤"
    #     gte <- "≥"
    #     paste(c(sprintf("Up Left: %s genes\nTP %s %s\nFDP %s %s", req(TP_FDP()$n2), gte,
    #                     req(TP_FDP()$TP2), lte, req(TP_FDP()$FDP2))))
    # })
    
    annotation <- reactive({ 
        A <- tibble::rownames_to_column(data.frame(df()), "VALUE")
        
        B = data()$annotation[c('affy_hg_u95av2','hgnc_symbol')]
        B
        
        B <- dplyr::left_join(A,B, by=c("VALUE"="affy_hg_u95av2"))[c('VALUE','hgnc_symbol')]
        rownames(B) <- B[['VALUE']]
        B <- B['hgnc_symbol']
        return(B)
    })
    
    output$volcanoplot <- renderPlotly({
        
        lte <- "≤"
        gte <- "≥"
        f <- list(
            size = 14,
            color = "#000000"
        )
        col = c("#33333333", "#FF0000", "#FF666633")
        cols <- rep(col[1], length(df()$fc))
        cols[c(selectedGenes()$sel1, selectedGenes()$sel2)] <- col[2]
        
        withProgress( message = "Please wait", { 
            
            plot_ly(data.frame(x = df()$fc, y=df()$logp), 
                    x = ~x, y = ~y, 
                    marker = list(size = 2,
                                  showlegend = FALSE,
                                  color = 'grey'), 
                    name = 'unselected',
                    type='scatter', source='A')%>%
                add_markers(x = selected_points()$x, y = selected_points()$y,
                            marker = list(
                                color = "red",
                                size = 6
                            ), 
                            name="selected") %>%
                layout(
                    showlegend = FALSE,
                    xaxis = list(title = "Fold change (log scale)", titlefont = f),
                    yaxis = list(title = "p-value (-log[10] scale)", titlefont = f), 
                    title = "Volcano plot",
                    # title = paste(c(sprintf("%s genes selected\nAt least %s true positives (FDP %s %s)",
                    #                         req(TP_FDP()$n12), req(TP_FDP()$TP12), lte, req(TP_FDP()$FDP12)))),
                    shapes = list(
                        list(
                            type = "line",
                            line = list(color = "orange", dash = "dot"),
                            x0 = xint(),
                            x1 = xint(),
                            y0 = 0,
                            y1 = 1,
                            yref = "paper"
                        )
                        ,
                        list(
                            type = "line",
                            line = list(color = "orange", dash = "dot"),
                            x0 = 0,
                            x1 = 1,
                            y0 = yint(),
                            y1 = yint(),
                            xref = "paper"
                        )
                        ,
                        list(
                            type = "line",
                            line = list(color = "orange", dash = "dot"),
                            x0 = xint2(),
                            x1 = xint2(),
                            y0=0,
                            y1 = 1,
                            yref = "paper"
                        )
                        # , 
                        # list(type = "rect",
                        #      fillcolor = "red", line = list(color = "red"), opacity = 0.2,
                        #      x0 = min(df()$fc)-0.5, x1 = xint2(), xref = "x",
                        #      y0 = yint(), y1 = max(df()$logp)+0.5 , yref = "y"), 
                        # list(type = "rect",
                        #      fillcolor = "red", line = list(color = "red"), opacity = 0.2,
                        #      x0 = xint(), x1 = max(df()$fc)+0.5, xref = "x",
                        #      y0 = yint(), y1 = max(df()$logp)+0.5 , yref = "y")
                    )
                ) %>% add_annotations(
                    x= 0,
                    y= 1,
                    xref = "paper",
                    yref = "paper",
                    text = paste(c(sprintf("<b>%s genes\nTP %s %s\nFDP %s %s</b>", req(TP_FDP()$n2), gte,
                                           req(TP_FDP()$TP2), lte, req(TP_FDP()$FDP2)))),
                    showarrow = F
                )%>% add_annotations(
                    x= 1,
                    y= 1,
                    xref = "paper",
                    yref = "paper",
                    text = paste(c(sprintf("<b>%s genes\nTP %s %s\nFDP %s %s</b>", req(TP_FDP()$n1), gte,
                                           req(TP_FDP()$TP1), lte, req(TP_FDP()$FDP1)))),
                    showarrow = F
                )%>% add_markers(
                    showlegend = FALSE,
                    text = annotation()[['hgnc_symbol']],
                    customdata = paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", annotation()[['hgnc_symbol']]))%>%
                onRender("
                    function(el) {
                        el.on('plotly_click', function(d) {
                            var url = d.points[0].customdata;
                            window.open(url);
                        });
                    }
                ") %>%
                
                
                config(editable = TRUE)
        })
    })
    
})
