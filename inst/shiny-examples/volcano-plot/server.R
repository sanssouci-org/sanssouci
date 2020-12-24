shinyServer(function(input, output, session) {
  
  source("function.R")
  
  options(shiny.maxRequestSize=1024^3)
  output$help <- renderUI({
    a("IIDEA help page", href = "https://pneuvial.github.io/sanssouci/articles/IIDEA.html", target= "_blank")
  })

  exampleData <- reactive({
    data <- list()
    data$matrix <- expr_ALL
    
    categ <- colnames(data$matrix)
    data$categ <- rep(1, length(categ))
    data$categ[which(categ == "NEG")] <- 0
    
    colnames(data$matrix) <- data$categ
    
    data$geneNames <- rownames(data$matrix)
    
    # data$biologicalFunc <- defaultBiologicalFunc(expr_ALL, expr_ALL_annotation)
    bioFun <- expr_ALL_GO
    stopifnot(nrow(bioFun) == nrow(data$matrix))  ## sanity check: dimensions
    ## make sure the ordering of probes (genes) 
    ## is the same for biological functions and expression data:
    mm <- match(rownames(bioFun), rownames(data$matrix))
    stopifnot(!any(is.na(mm)))
    data$biologicalFunc <- bioFun[mm, ]
    rm(bioFun)
    return(data)
  })
  
  output$downloadExampleData <- downloadHandler(
    filename = function() {
      paste("ExampleData", "zip", sep=".")
    },
    content = function(fname) {
      fs <- c()
      tmpdir <- tempdir()
      setwd(tempdir())
      paths <- c("expressData.csv", "biologicalFunction.csv")
      write.csv(exampleData()$matrix, paths[1])
      write.csv(exampleData()$biologicalFunc, paths[2])
      zip(zipfile=fname, files=paths)
    },
    contentType = "application/zip"
  )
  
  
  data <- 
    reactive (
      # eventReactive(input$buttonData,
      {
        # output$CheckData <- renderUI({
        #   tags$span(style="color:blue", "Please, upload your data")
        # })
        data <- list()
        if (input$checkboxDemo){
          data$matrix <- expr_ALL
          
          categ <- colnames(data$matrix)
          data$categ <- rep(1, length(categ))
          data$categ[which(categ == "NEG")] <- 0
          
          data$geneNames <- rownames(data$matrix)
          
          # data$biologicalFunc <- defaultBiologicalFunc(expr_ALL, expr_ALL_annotation)
          bioFun <- expr_ALL_GO
          stopifnot(nrow(bioFun) == nrow(data$matrix))  ## sanity check: dimensions
          ## make sure the ordering of probes (genes) 
          ## is the same for biological functions and expression data:
          mm <- match(rownames(bioFun), rownames(data$matrix))
          stopifnot(!any(is.na(mm)))
          data$biologicalFunc <- bioFun[mm, ]
          rm(bioFun)
          # output$CheckData <- renderUI({
          #   actionButton("buttonValidate", "Perform calibration!")
          # })
          
        } else {
          req(input$fileData)
          file <- req(input$fileData)
          data$matrix <- read.csv(file = file$datapath, row.names = 1, check.names=FALSE)
          
          colnam <- colnames(data$matrix) #Control colnames(matrix)
          valueColnam <- sort(unique(colnam))
          if(length(valueColnam)==2){
            if(valueColnam != c("0","1")){
              colnam[ colnam == valueColnam[1] ] <- 0
              colnam[ colnam == valueColnam[2] ] <- 1
              colnames(data$matrix) <- colnam
              output$errorInput <- renderUI({
                tags$span(style="color:orange", paste("The colnames of your matrix does not contains 0 or 1.We consider that", valueColnam[1] ,"becomes 0 and ", valueColnam[2]," becomes 1"))
              })
            }
          }else{
            output$errorInput <- renderUI({
              tags$span(style="color:red", "The column names of your data contains more (or less) than 2 categories. Please use {0, 1} for the colnames of your matrix.")
            })
          }
          data$categ <- colnames(data$matrix)
          
          data$geneNames <- rownames(data$matrix)
          
          req(input$fileGroup)
          fileGroup <- req(input$fileGroup)
          biofun <- read.csv(file = fileGroup$datapath, row.names = 1, check.names=FALSE)
          stopifnot(nrow(bioFun) == nrow(data$matrix))  ## sanity check: dimensions
          ## make sure the ordering of probes (genes) 
          ## is the same for biological functions and expression data:
          mm <- match(rownames(bioFun), rownames(data$matrix))
          stopifnot(!any(is.na(mm)))
          data$biologicalFunc <- bioFun[mm, ]
          rm(bioFun)
          
          # output$CheckData <- renderUI({
          #   actionButton("buttonValidate", "Perform calibration!")
          # })
          
        }
        return(data)
      }
    )
  output$inputK <- renderUI({
    numericInput("valueK", 
                 label = "K (size of reference family)", 
                 value = ifelse(input$refFamily == "Beta", 
                                round(2*nrow(data()$matrix)/100),
                                nrow(data()$matrix)),
                 min = 1,
                 max = nrow(data()$matrix))
  })
  
  
  alpha <- eventReactive(input$buttonValidate, {
    req(1 - input$sliderConfLevel/100)
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
                 refFamily = refFamily(), alternative = alternative(), 
                 K = numK()
    )
  })
  
  df <- reactive({
    dex <- rowWelchTests(data()$matrix, data()$categ)
    pval <- dex[["p.value"]]
    logp <- -log10(pval)
    fc <- dex$meanDiff
    adjp <- p.adjust(pval, method = "BH")
    return(list(logp = logp, fc = fc, adjp = adjp, pval = pval))
  })
  
  
  vertical <- reactive({event_data("plotly_relayout", source='A')})
  
  xint <- reactiveVal(0.5)
  observeEvent(vertical()[["shapes[0].x0"]], {
    if(input$symetric){
      newValue <- vertical()[["shapes[0].x0"]]
      xint(newValue)
      xint2(-newValue)
    } else {
      newValue <- vertical()[["shapes[0].x0"]]   
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
                    round(xint(), digits = 3), 
                    round(xint2(), digits = 3),  
                    formatC(10^(-yint()), format = "e", digits = 3)
    )))
  })
  
  selectedGenes <- reactive({
    req(df())
    ## gene selections
    sel1 <- which(df()$logp >= yint() & df()$fc >= xint()) 
    sel2 <- which(df()$logp >= yint() & df()$fc <= xint2())
    sel12 <- sort(union(sel1,sel2))
    return(list(sel1 = sel1, sel2 = sel2, sel12 = sel12))
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
    
    return(list(n1 = n1, TP1 = TP1, FDP1 = FDP1, n2 = n2, TP2 = TP2, FDP2 = FDP2, n12 = n12, TP12 = TP12, FDP12 = FDP12))
    
  })
  
  selected_points <- reactive({
    list(x = df()$fc[selectedGenes()$sel12], 
         y = df()$logp[selectedGenes()$sel12])
  })
  
  tableResult <- reactiveVal(data.frame(
    Selection = c("Threshold selection"))) #Initialization
  
  baseTable <- reactive({
    data.frame(`Selection` = c("Threshold selection"), 
               "# genes" = c(TP_FDP()$n12),
               "TP≥" = as.integer(c(TP_FDP()$TP12)), 
               "FDP≤" = c(TP_FDP()$FDP12),
               check.names = FALSE)
  })
  
  observeEvent(TP_FDP(),{ # When threshold change
    
    bottomTable <- tableResult() %>% 
      filter(Selection != "Upper right") %>% 
      filter(Selection != "Threshold selection") %>% 
      filter(Selection != "Upper left")
    upperTable <- baseTable()
    newValue <- rbind(upperTable, bottomTable)
    tableResult(newValue)
  })
  observeEvent(d(),{  # When user selects a new group of points
    req(calcBoundSelection())
    vectorGene <- names(df()$pval[manuelSelected()])
    url <- UrlStringdbGrah(vectorGene)
    n <- dim(tableResult())[1]
    newValue <- rbind(tableResult(), c(paste('<a target="_blanck" href="', url, '" >User selection ',n, '</a>', sep=""), 
                                       calcBoundSelection()$n, 
                                       calcBoundSelection()$TP,
                                       calcBoundSelection()$FDP))
    tableResult(newValue)
  })
  observeEvent(input$resetCSV, { # to clean printed table
    newValue <- baseTable()
    tableResult(newValue)
  })
  observeEvent(data(), { # to clean printed table
    newValue <- baseTable()
    tableResult(newValue)
  })
  
  
  
  textQtableBounds <- reactive({
    paste("This table prints your post-hoc bounds for your selections.", 
          "FDP : False discovery Proportion.", 
          "TP : True positive", 
          "For example, the selection called Upper Left have at leat", 
          tableResult()[1, "TP≥"],
          " true positives and max ",
          tableResult()[1, "FDP≤"],
          " False discovery proportion.")})
  
  output$OutQtableBounds <- renderUI({
    tab <- tableResult()
    msg <- "The table below prints post hoc bounds for user selections. "
    if (nrow(tab)>0) {
      msg <- paste(msg, "For example, the selection called",
                       tab[1, 1], 
                       "contains at least", 
                       tab[1, "TP≥"],
                       " true positives (TP) and its False Discovery Proportion (FDP) is less than",
                       tab[1, "FDP≤"])
    }
    popify(el = bsButton("QtableBounds", label = "", icon = icon("question"), style = "info", size = "extra-small"), 
           title = "Data", 
           content = msg,
           trigger='hover')
  })
  
  # addPopover(session, "QtableBounds", "Data", content = textQtableBounds(), trigger = 'focus')
  
  output$tableBounds <- renderDT({
    req(TP_FDP())
    
    tableResult()
    
  }, selection = list(mode = 'single', selectable = -(1)) , escape = FALSE )
  
  
  lineAdjp <- reactive({ # value for 
    listLog <- c()
    for (i in c(0.5,0.25,0.1,0.05,0.025,0.01,0.001,0.0001)){
      min05 <- names(which.min(abs(df()$adjp-i)))
      y05 <- df()$logp[[min05]]
      listLog <- c(listLog, y05)
    }
    return(listLog)
  })
  
  output$lineAdjp <- renderPrint({list(lineAdjp(), class(lineAdjp()))})
  
  thr_yaxis <- reactive({
    thrYaxis(thr = cal()$thr, maxlogp=max(df()$logp))
  })
  
  yaxis <- reactive({
    f <- list(
      size = 14,
      color = "#000000"
    )
    
    yaxis <- switch(input$choiceYaxis, 
                    "pval" = list(
                      title = "p-value (-log[10] scale)", 
                      titlefont = f
                    ),
                    "adjPval" = list(
                      title = "Adjusted p-value (-log[10] scale)", 
                      titlefont = f,
                      autotick = FALSE,
                      tickmode = "array",
                      tickvals = lineAdjp(),
                      ticktext = c(0.5,0.25,0.1,0.05,0.025,0.01,0.001,0.0001)
                    ),
                    "thr" = list(
                      title = "Calibration thresholds (-log[10] scale)", 
                      titlefont = f, 
                      autotick = FALSE,
                      tickmode = "array",
                      tickvals = thr_yaxis()$pvalue,
                      ticktext = thr_yaxis()$num
                    )
    )
    return(yaxis)
  })
  
  thrLine <- reactive({
    list(
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
    )
  })
  
  
  
  posteriori <- eventReactive(input$buttonValidate, {
    f <- list(
      size = 14,
      color = "#000000"
    )
    lte <- "≤"
    gte <- "≥"
    plot_ly(data.frame(x = df()$fc, y=df()$logp), 
            x = ~x, y = ~y, 
            marker = list(size = 2,
                          # showlegend = FALSE,
                          color = 'grey'), 
            name = 'unselected',
            type='scattergl', mode = "markers", source='A'
            # , text = annotation()[['nameGene']],
            # customdata = paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", annotation()[['nameGene']])
            # , height = 600
    ) %>%
      add_markers(x = selected_points()$x, y = selected_points()$y,
                  marker = list(
                    color = "red",
                    size = 6
                  ),
                  name ="selected") %>%
      layout(
        # showlegend = FALSE,
        xaxis = list(title = "Fold change (log scale)", titlefont = f),
        yaxis = isolate(yaxis()),
        title = "",
        shapes = isolate(thrLine()),
        dragmode = "select" ) %>%
      # add_annotations(
      #     x = 0,
      #     y = 1,
      #     xref = "paper",
      #     yref = "paper",
      #     text = paste(c(sprintf(
      #         "<b>%s genes\nTP %s %s\nFDP %s %s</b>", 
      #         req(TP_FDP()$n2), gte,
      #         req(TP_FDP()$TP2), lte, 
      #         req(TP_FDP()$FDP2)))),
      #     showarrow = F
    # ) %>% add_annotations(
    #     x = 1,
    #     y = 1,
    #     xref = "paper",
    #     yref = "paper",
    #     text = paste(c(sprintf(
    #         "<b>%s genes\nTP %s %s\nFDP %s %s</b>", 
    #         req(TP_FDP()$n1), gte,
    #         req(TP_FDP()$TP1), lte, 
    #         req(TP_FDP()$FDP1)))),
    #     showarrow = FALSE
    # ) %>% 
    onRender("
                  function(el) {
                      el.on('plotly_click', function(d) {
                          var url = d.points[0].customdata;
                          window.open(url);
                      });
                  }
        ") %>%
      event_register("plotly_selecting") %>%
      config(editable = TRUE)%>%
      toWebGL()
  })
  
  output$volcanoplotPosteriori <- renderPlotly({
    withProgress( message = "Please wait", {
      posteriori()
    })
  })
  
  
  
  observeEvent( vertical(), { #when threshold moved and red selected point change
    
    plotlyProxy("volcanoplotPosteriori", session) %>%
      plotlyProxyInvoke("deleteTraces",1)%>%
      plotlyProxyInvoke(
        "addTraces",
        list(
          x = unname(selected_points()$x),
          y = unname(selected_points()$y),
          type = "scatter",
          mode = "markers",
          line = list(color = "red"), 
          name = "selected", 
          showlegend = TRUE
        ), 1
      )
    
    
  })
  
  
  observeEvent(input$choiceYaxis, { #when we choose a different y axis 
    plotlyProxy("volcanoplotPosteriori", session) %>%
      plotlyProxyInvoke("relayout", list(yaxis = yaxis()))
    
    
  })
  
  observeEvent(vertical(), {  # for moving of orange threshold
    plotlyProxy("volcanoplotPosteriori", session)%>%
      plotlyProxyInvoke("relayout", 
                        list(shapes = thrLine()))
  })
  
  
  
  d <- reactive({ event_data("plotly_selected", source='A')})
  
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
  
  
  calcBoundSelection <- reactive({ #calculate bounds of selected genes 
    req(manuelSelected())
    calcBounds(df()$pval[manuelSelected()], thr = cal()$thr)
  })
  
  output$downloadData <- downloadHandler( #download csv of user selection
    filename = function() {
      paste("SelectionList", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(tableCSV(), file)
    }
  )
  
  output$curveMaxFPBoth <- renderPlotly({
    plotMaxFP(pval = df()$pval[selectedGenes()$sel12], thr = cal()$thr) + 
      ggtitle("Upper Left + right")
  })
  
  output$curveMaxFPSelect <- renderPlotly({
    if(length(userDTselectPost()) == 1){
      plotMaxFP(pval = df()$pval[selectionUserRe()$sel], thr = cal()$thr) + 
        ggtitle(userDTselectPost())
    }else{
      plotMaxFP(pval = df()$pval[manuelSelected()], thr = cal()$thr) + 
        ggtitle("User selection")
    }
  })
  
  
  
  
  userDTselectPost <- reactive({
    href <- tableResult()[input$tableBounds_rows_selected, "Selection"]
    str_remove_all(str_remove_all(href, "<a(.*?)>"), "(</a>)")
  })
  
  selectionUserRe <- reactive({
    vect <- tableCSV()[,userDTselectPost()]
    sel <- which(vect == 1)
    list(sel = sel)
  })
  
  observeEvent(userDTselectPost(), {
    if(length(userDTselectPost()) == 1){
      plotlyProxy("volcanoplotPosteriori", session) %>%
        plotlyProxyInvoke("deleteTraces", 2)
      plotlyProxy("volcanoplotPosteriori", session) %>%
        plotlyProxyInvoke(
          "addTraces",
          list(
            x = unname(df()$fc[selectionUserRe()$sel]),
            y = unname(df()$logp[selectionUserRe()$sel]),
            type = "scattergl",
            mode = "markers",
            line = list(color = "blue"),
            name = userDTselectPost()
          ), 2
        )
    } else {
      plotlyProxy("volcanoplotPosteriori", session) %>%
        plotlyProxyInvoke("deleteTraces",2)
    }
  })
  
  
  
  
  
  ## biological function 
  tableBoundsGroup <- reactive({
    req(data())
    
    boundGroup(df(), 
               data()$biologicalFunc, 
               thr = cal()$thr,
               nameFunctions = colnames(data()$biologicalFunc))
  })
  
  output$OutQtableBoundsGroup <- renderUI({
    popify(el = bsButton("QtableBoundsGroup", label = "", icon = icon("question"), style = "info", size = "extra-small"), 
           title = "Data", content = paste("This table prints your post-hoc bounds for your selections."  ,
                                           "FDP : False discovery Proportion.",
                                           "TP : True positive", 
                                           "For example, the selection called Upper Left have at leat", 
                                           tableBoundsGroup()[1, "TP≥"],
                                           " true positives and max ",
                                           tableBoundsGroup()[1, "FDP≤"]*100,
                                           "% False discovery proportion.")
           , trigger='focus')
  })
  
  output$tableBoundsGroup <- renderDT({
    tableBoundsGroup()
  }, selection = 'single')
  
  # output$choiceGroupUI <- renderUI({
  #   selectInput("choiceGroup", label = "Gene set", 
  #               choices = c("Select a gene set", colnames(data()$biologicalFunc))
  #   )
  # })
  
  userDTselectPrio <- reactive({
    tableBoundsGroup()[input$tableBoundsGroup_rows_selected,"Name"]
  })
  
  selectionGroup <- reactive({
    req(data())
    req(df())
    
    group <- req(userDTselectPrio())
    bioFun <- data()$biologicalFunc
    ids <- which(bioFun[, group] == 1)
    list(sel = ids)
  })
  
  priori <- eventReactive(input$buttonValidate, {
    f <- list(
      size = 14,
      color = "#000000"
    )
    lte <- "≤"
    gte <- "≥"
    plot_ly(data.frame(x = df()$fc, y=df()$logp), 
            x = ~x, y = ~y, 
            marker = list(size = 2,
                          # showlegend = FALSE,
                          color = 'grey'), 
            name = 'genes',
            type='scattergl', mode = "markers", source='B'
            ,
            text = data()$geneName,
            customdata = paste0("http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=", data()$geneName)
            # , height = 600
    )%>% 
      layout(
        showlegend = TRUE,
        xaxis = list(title = "Fold change (log scale)", titlefont = f),
        yaxis = isolate(yaxis()),
        title = "",
        # shapes = isolate(thrLine()),
        dragmode = "select" )%>%
      onRender("
                  function(el) {
                      el.on('plotly_click', function(d) {
                          var url = d.points[0].customdata;
                          window.open(url);
                      });
                  }
              ") %>%
      event_register("plotly_selecting") %>%
      config(editable = TRUE) %>%
      toWebGL() 
    
  })
  
  
  output$volcanoplotPriori <- renderPlotly({
    withProgress( message = "Please wait", {
      priori()
    })
  })
  
  observeEvent(input$choiceYaxis, { #when we choose a different y axis 
    plotlyProxy("volcanoplotPriori", session) %>%
      plotlyProxyInvoke("relayout", list(yaxis = yaxis()))
    
    
  })
  
  
  
  observeEvent(userDTselectPrio(), {
    if(length(userDTselectPrio()) == 1){
      plotlyProxy("volcanoplotPriori", session) %>%
        plotlyProxyInvoke("deleteTraces", 1)
      plotlyProxy("volcanoplotPriori", session) %>%
        plotlyProxyInvoke(
          "addTraces",
          list(
            x = unname(df()$fc[selectionGroup()$sel]),
            y = unname(df()$logp[selectionGroup()$sel]),
            type = "scattergl",
            mode = "markers",
            line = list(color = "blue"),
            name = userDTselectPrio()
          )
        )
    } else {
      plotlyProxy("volcanoplotPriori", session) %>%
        plotlyProxyInvoke("deleteTraces",1)
    }
  })
  
  output$curveMaxFPGroup <- renderPlotly({
    req(selectionGroup())
    plotMaxFP(pval = df()$pval[selectionGroup()$sel], thr = cal()$thr) + 
      ggtitle(userDTselectPrio()) 
  })
  
  
  
})
