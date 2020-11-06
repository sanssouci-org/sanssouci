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
    library(dplyr)
    
    source("function.R")
    
    options(shiny.maxRequestSize=1024^3)
    
    data <- 
        reactive (
            # eventReactive(input$buttonData,
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
    output$inputK <- renderUI({
        numericInput("valueK", label="K - il est calculé qu'une fois que les données sont validées", value = nrow(data()$matrix), max=nrow(data()$matrix))
    })
    
    output$datatable <- renderTable({head(data()$matrix)})
    
    alpha <- eventReactive(input$buttonValidate, {
        req(input$sliderAlpha)
    })
    numB <-eventReactive(input$buttonValidate, {
        req(input$numB)
    })
    refFamily <- eventReactive(input$buttonValidate,{
        req(input$refFamily)
    })
    alternative <- eventReactive(input$buttonValidate, {
        req(input$alternative)
    })
    numK <- eventReactive(input$buttonValidate, {
        req(input$valueK)
    })
    cal <- reactive({
        calibrateJER(data()$matrix, categ = data()$categ, 
                     B = numB(), alpha = alpha(), 
                     refFamily=refFamily(), alternative = alternative(), 
                     K = numK()
        )
    })
    
    df <- reactive({
        dex <- rowWelchTests(data()$matrix, data()$categ)
        pval <- dex[["p.value"]]
        logp <- -log10(pval)
        fc <- dex$meanDiff
        adjp <- p.adjust(pval, method = "BH")
        return(list(logp=logp, fc=fc, adjp=adjp, pval=pval))
    })
    
    
    
    
    
    
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
        paste(c(sprintf("Foldchange: \n right: %s \nleft: %s \n P_values: %s", 
                        round(xint(), digits = 3), round(xint2(), digits = 3),  formatC(10^(-yint()), format = "e", digits = 3)
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
    
    tableResult <- reactiveVal( data.frame(Selection=c("Upper Right", "Upper Left", "Both parts"))) #Initialization
    observeEvent(TP_FDP(),{ #When threshold change
        
        bottomTable <- tableResult() %>% 
            filter(Selection != "Upper Right") %>% 
            filter(Selection != "Both parts") %>% 
            filter(Selection!="Upper Left")
        upperTable <- data.frame(`Selection`=c("Upper Right", "Upper Left", "Both parts"), 
                               "Nb of genes"=c(TP_FDP()$n1, TP_FDP()$n2, TP_FDP()$n12),
                               "TP≥"=as.integer(c(TP_FDP()$TP1, TP_FDP()$TP2, TP_FDP()$TP12)), 
                               "FDP≤" = c(TP_FDP()$FDP1, TP_FDP()$FDP2, TP_FDP()$FDP12))
        newValue <- rbind(upperTable, bottomTable)
        tableResult(newValue)
    })
    observeEvent(d(),{ #When user select a new group of points
        req(calcBoudSelection())
        n <- dim(tableResult())[1]-2
        newValue <- rbind(tableResult(), c(paste("User selection",n), 
                                           calcBoudSelection()$n, 
                                           calcBoudSelection()$TP,
                                           calcBoudSelection()$FDP))
        tableResult(newValue)
    })
    observeEvent(input$resetCSV, { # to clean printed table
        newValue <- data.frame(`Selection`=c("Upper Right", "Upper Left", "Both parts"), 
                               "Nb of genes"=c(TP_FDP()$n1, TP_FDP()$n2, TP_FDP()$n12),
                               "TP≥"=as.integer(c(TP_FDP()$TP1, TP_FDP()$TP2, TP_FDP()$TP12)), 
                               "FDP≤" = c(TP_FDP()$FDP1, TP_FDP()$FDP2, TP_FDP()$FDP12))
        tableResult(newValue)
    })
    
    
    
    output$tableBounds <- renderTable({
        req(TP_FDP())
        
        tableResult()
        
    })
    
    
    lineAdjp <- reactive({ # value for 
        listLog <- c()
        for (i in c(0.5,0.25,0.1,0.05,0.025,0.01,0.001,0.0001)){
            min05 <- names(which.min(abs(df()$adjp-i)))
            y05 <- df()$logp[min05]
            listLog <- c(listLog, y05)
            print(y05)
        }
        return(listLog)
    })

    
    annotation <- reactive({ 
        A <- tibble::rownames_to_column(data.frame(df()), "VALUE")
        
        B = data()$annotation[c('affy_hg_u95av2','hgnc_symbol')]
        B
        
        B <- dplyr::left_join(A,B, by=c("VALUE"="affy_hg_u95av2"))[c('VALUE','hgnc_symbol')]
        rownames(B) <- B[['VALUE']]
        B <- B['hgnc_symbol']
        return(B)
    })
    
    yaxis <- reactive({
        f <- list(
            size = 14,
            color = "#000000"
        )
        
        yaxis <- switch(input$choiceYaxis, 
                        "pval"= list(
                            title = "p-value (-log[10] scale)", 
                            titlefont = f
                        ),
                        "adjPval"=list(
                            title = "Adjusted p_value (-log[10] scale)", 
                            titlefont = f,
                            autotick = FALSE,
                            tickmode = "array",
                            tickvals = lineAdjp(),
                            ticktext = c(0.5,0.25,0.1,0.05,0.025,0.01,0.001,0.0001)
                        ),
                        "thr"=list(
                            title = "Calibration thr (-log[10] scale)", 
                            titlefont = f, 
                            autotick = FALSE,
                            tickmode = "array",
                            tickvals = -log10(cal()$thr)[1:10],
                            ticktext = as.character(1:10)
                        )
        )
        return(yaxis)
    })
    
    
    output$volcanoplot <- renderPlotly({
        
        f <- list(
            size = 14,
            color = "#000000"
        )
        lte <- "≤"
        gte <- "≥"
        
        
       
        col = c("#33333333", "#FF0000", "#FF666633")
        cols <- rep(col[1], length(df()$fc))
        cols[c(selectedGenes()$sel1, selectedGenes()$sel2)] <- col[2]
        
        withProgress( message = "Please wait", { 
            
            plot_ly(data.frame(x = df()$fc, y=df()$logp), 
                    x = ~x, y = ~y, 
                    marker = list(size = 2,
                                  # showlegend = FALSE,
                                  color = 'grey'), 
                    name = 'unselected',
                    type='scatter', mode = "markers", source='A'
                    # , height = 600
            )%>%
                add_markers(x = selected_points()$x, y = selected_points()$y,
                            marker = list(
                                color = "red",
                                size = 6
                            ),
                            name="selected") %>%
                # add_segments(x = min(df()$fc), xend = max(df()$fc), 
                #              y = lineAdjp()$logp0.10, yend = lineAdjp()$logp0.10, 
                #              name="Adjusted pvalues = 0.1", legendgroup = "Adjusted pvalues",
                #              opacity=0.3)%>%
                # add_segments(x = min(df()$fc), xend = max(df()$fc),
                #              y = lineAdjp()$logp0.05, yend = lineAdjp()$logp0.05, 
                #              name="Adjusted pvalues = 0.05", legendgroup = "Adjusted pvalues",
                #              opacity=0.3)%>%
                # add_segments(x = min(df()$fc), xend = max(df()$fc), 
                #              y = lineAdjp()$logp0.01, yend = lineAdjp()$logp0.01, 
                #              name="Adjusted pvalues = 0.01", legendgroup = "Adjusted pvalues",
                #              opacity=0.3)%>%
                layout(
                    # showlegend = FALSE,
                    xaxis = list(title = "Fold change (log scale)", titlefont = f),
                    yaxis = yaxis(), 
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
                        # list(
                        #     type = "line",
                        #     line = list(color = "grey"),
                        #     x0 = 0,
                        #     x1 = 1,
                        #     y0 = lineAdjp()$logp0.10,
                        #     y1 = lineAdjp()$logp0.10,
                        #     xref = "paper"
                        # ),
                        # list(
                        #     type = "line",
                        #     line = list(color = "grey"),
                        #     x0 = 0,
                        #     x1 = 1,
                        #     y0 = lineAdjp()$logp0.05,
                        #     y1 = lineAdjp()$logp0.05,
                        #     xref = "paper"
                        # ),
                        # list(
                        #     type = "line",
                        #     line = list(color = "grey"),
                        #     x0 = 0,
                        #     x1 = 1,
                        #     y0 = lineAdjp()$logp0.01,
                        #     y1 = lineAdjp()$logp0.01,
                        #     xref = "paper"
                        # )
                    ), 
                    dragmode = "select" ) %>%
                add_annotations(
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
                event_register("plotly_selecting") %>%
                config(editable = TRUE)
        })
    })
    
    d <- reactive({ event_data("plotly_selected")})
    
    manuelSelected <- reactive({
        req(d())
        d()[['pointNumber']][which(d()[['curveNumber']]==0)]+1
    })
    
    tableCSV <- reactiveVal(data.frame(row.names = rownames(expr_ALL)))
    observeEvent(input$resetCSV,{
        tableCSV(data.frame(row.names = rownames(expr_ALL)))
    })
    observeEvent(d(), { 
        vecteur <- rep(0, dim(data()$matrix)[1])
        vecteur[manuelSelected()] <- 1
        nameCol <- colnames(tableCSV())
        df <- cbind(tableCSV(), selection2=vecteur)
        colnames(df) <- c(nameCol, paste("User selection", length(nameCol)+1))
        tableCSV(df)
    })
    
    
    calcBoudSelection <- reactive({ #calculate bounds of selected genes 
        req(manuelSelected())
        calcBounds(df()$pval, selection=manuelSelected(), thr=cal()$thr)
    })
    
    output$downloadData <- downloadHandler( #download csv of user selection
        filename = function() {
            paste("SelectionList", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
            write.csv(tableCSV(), file)
        }
    )
    
})
