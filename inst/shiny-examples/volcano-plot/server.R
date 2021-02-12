shinyServer(function(input, output, session) {
  
  source("function.R")
  
  options(shiny.maxRequestSize=1024^3)
  output$help <- renderUI({
    a("IIDEA help page", href = "https://pneuvial.github.io/sanssouci/articles/IIDEA.html", target= "_blank")
  })
  
  geo2kegg <- reactive({
    withProgress(message = "Load GSEABenchmarkeR data set ... ", {
      t1 <- Sys.time()
      data <- R.cache::memoizedCall(GSEABenchmarkeR::loadEData,"geo2kegg")
      t2 <- Sys.time()
      print(paste("Load GSEABenchmarkeR data set:",difftime(t2, t1)))
      setProgress(value = 1, detail = "Done")
      return(data)
    })
  })
  
  nameGeo2kegg <- reactive({namedGeo2kegg(geo2kegg())})
  
  output$choiceGSEAUI <- renderUI({
    selectInput("choiceGSEA", label = "Choose a gene data set", 
                choices = c('sansSouci example data with gene set'='OurData', nameGeo2kegg()))
  })
  
  
  
  exampleData <- reactive({
    withProgress(value = 0, message = "Preparation for downloaded data ...", {
      
      setProgress(value = 0.3, detail ="Download data ... ")
      data <- list()
      data$matrix <- expr_ALL
      
      setProgress(value = 0.5, detail = "Cleaning data ... ")
      
      categ <- colnames(data$matrix)
      data$categ <- rep(1, length(categ))
      data$categ[which(categ == "NEG")] <- 0
      
      colnames(data$matrix) <- data$categ
      
      data$geneNames <- base::rownames(data$matrix)
      
      setProgress(value = 0.6, detail = "Preparation of gene set data ...  ")
      # data$biologicalFunc <- defaultBiologicalFunc(expr_ALL, expr_ALL_annotation)
      bioFun <- expr_ALL_GO
      stopifnot(nrow(bioFun) == nrow(data$matrix))  ## sanity check: dimensions
      ## make sure the ordering of probes (genes) 
      ## is the same for biological functions and expression data:
      mm <- match(base::rownames(bioFun), base::rownames(data$matrix))
      stopifnot(!any(is.na(mm)))
      data$biologicalFunc <- bioFun[mm, ]
      rm(bioFun)
      setProgress(value = 1, detail = "Done")
      return(data)
    })
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
  
  fileData <- reactiveVal(NULL)
  observe({
    newValue <- req(input$fileData)
    fileData(newValue)
  })
  
  fileGroup <- reactiveVal(NULL)
  observe({
    newValue <- input$fileGroup
    fileGroup(newValue)
  })
  
  data <- 
    # reactive (
    eventReactive(input$buttonValidate,
                  {
                    withProgress(value = 0, message = "Upload Data... ", {
                      data <- list()
                      if (input$checkboxDemo){
                        if(req(input$choiceGSEA)=='OurData'){
                          setProgress(value = 0.4, detail = "SansSoucis data set ...")
                          data$matrix <- expr_ALL
                          
                          categ <- colnames(data$matrix)
                          data$categ <- rep(1, length(categ))
                          data$categ[which(categ == "NEG")] <- 0
                          
                          data$geneNames <- rownames(data$matrix)
                          
                          setProgress(value = 0.7, detail = "Preparation of gene set data ...  ")
                          
                          # data$biologicalFunc <- defaultBiologicalFunc(expr_ALL, expr_ALL_annotation)
                          bioFun <- expr_ALL_GO
                          stopifnot(nrow(bioFun) == nrow(data$matrix))  ## sanity check: dimensions
                          ## make sure the ordering of probes (genes) 
                          ## is the same for biological functions and expression data:
                          
                          mm <- match(base::rownames(bioFun), base::rownames(data$matrix))
                          stopifnot(!any(is.na(mm)))
                          data$biologicalFunc <- bioFun[mm, ]
                          rm(bioFun)
                          # output$CheckData <- renderUI({
                          #   actionButton("buttonValidate", "Perform calibration!")
                          # })
                          data$boolValidation <- TRUE 
                        } else {
                          setProgress(value = 0.4, detail = "GSEA data set ...")
                          rawData <- R.cache::memoizedCall(maPreproc,geo2kegg()[input$choiceGSEA])[[1]]
                          # rawData <- maPreproc(geo2kegg()[input$choiceGSEA])[[1]]
                          
                          data$matrix <- SummarizedExperiment::assays(rawData)$exprs
                          
                          cats <- SummarizedExperiment::colData(rawData)
                          ww <- match(cats$Sample, base::colnames(data$matrix))
                          categ <- cats$GROUP[ww]
                          # colnames(data$matrix) <- categ
                          data$categ <- categ
                          setProgress(value = 0.7, detail = "GSEA data set ...")
                          
                          data$boolValidation <- TRUE
                          
                          print(input$choiceGSEA)
                        }
                        
                      } else { # if user use his own data
                        req(input$fileData)
                        req(fileData())
                        file <- req(fileData())
                        setProgress(value = 0.1, detail = "Read csv ...")
                        data$matrix <- read.csv(file = file$datapath, row.names = 1, check.names=FALSE)
                        
                        data$boolDegrade <- (length(data$matrix) == 2 & all(sort(colnames(data$matrix)) == sort(c("fc","p.value"))))
                        setProgress(value = 0.2, detail = "Test data ...")
                        if (data$boolDegrade){
                          setProgress(value = 0.4)
                          data$df <- data$matrix
                        } else {
                          setProgress(value = 0.4, detail = "Clean data set ...")
                          clean <- cleanMatrix(data$matrix)
                          
                          if(clean$boolValidation){
                            data$matrix = clean$data
                          }else{
                            data$matrix=NULL
                          }
                          
                          data$matrix.color <- clean$color
                          data$matrix.text <- clean$text
                          
                          data$categ <- colnames(data$matrix)
                          
                          
                          data$geneNames <- base::rownames(data$matrix)
                          data$boolValidation <- clean$boolValidation
                        }
                        
                        setProgress(value = 0.7, detail = "Preparation of gene set data ...")
                        # req(input$fileGroup)
                        fileGroup <- fileGroup()
                        if (!is.null(fileGroup)){
                          setProgress(value = 0.75, detail = "Read gene set data ...")
                          bioFun <- read.csv(file = fileGroup$datapath, row.names = 1, check.names=FALSE)
                          setProgress(value = 0.8, detail = "Cleaning of gene set data ...")
                          cleanBio <- cleanBiofun(bioFun)
                          
                          data$bioFun.color <- cleanBio$color
                          data$bioFun.text <- cleanBio$text
                          
                          
                          
                          
                          bioFun <- cleanBio$biofun
                          setProgress(value = 0.9, detail = "Matching ...")
                          
                          matchBio <- matchMatrixBiofun(matrixFunc = data$matrix, biofun = bioFun)
                          if(matchBio$boolValidation & cleanBio$boolValidation){
                            data$biologicalFunc <- matchBio$biofun
                          } else {
                            data$biologicalFunc <- NULL
                          }
                          
                          data$match.color <- matchBio$color
                          data$match.text <- matchBio$text
                          
                          
                          
                          
                          rm(bioFun)
                        }
                        
                        if (data$boolDegrade){
                          data$matrix <- NULL
                        } 
                        
                        
                        
                      }
                      setProgress(value = 1, detail = "Done")
                      return(data)
                    })
                  }
)

output$errorInput <- renderUI({
  tags$span(style= req(data()$matrix.color), paste(req(data()$matrix.text)))
})

output$errorBioMatrix <- renderUI({
  tags$span(style=req(data()$bioFun.color), paste(req(data()$bioFun.text)))
})
output$errorMatch <- renderUI({
  tags$span(style=req(data()$match.color), paste(req(data()$match.text)))
})

observeEvent(input$resetInputData, {
  fileData(NULL)
  reset("fileData")
})

observeEvent(input$resetInputGroup, {
  fileGroup(NULL)
  reset("fileGroup")
})


output$inputK <- renderUI({
  req(data()$matrix)
  numericInput("valueK", 
               label = "K (size of reference family)", 
               value = ifelse(input$refFamily == "Beta", 
                              round(2*nrow(data()$matrix)/100),
                              nrow(data()$matrix)),
               min = 1,
               max = nrow(data()$matrix))
})


alpha <- reactiveVal(0.1) #Initialization
observeEvent(input$buttonValidate,{ # When Run is clicked
  
  newValue <- req(1 - input$sliderConfLevel/100)
  alpha(newValue)
})
# alpha <-eventReactive(input$buttonValidate, {
#   req(1 - input$sliderConfLevel/100)
# })
numB <- reactiveVal(1000) #Initialization
observeEvent(input$buttonValidate, {
  newValue <- req(input$numB)
  numB(newValue)
})
# numB <-eventReactive(input$buttonValidate, {
#   req(input$numB)
# })
refFamily <- reactiveVal("Simes")
observeEvent(input$buttonValidate, {
  newValue <- req(input$refFamily)
  refFamily(newValue)
})
# refFamily <- eventReactive(input$buttonValidate, {
#   req(input$refFamily)
# })
alternative <- reactiveVal("two.sided")
observeEvent(input$buttonValidate, {
  newValue <- req(input$alternative)
  alternative(newValue)
})
# alternative <-eventReactive(input$buttonValidate, {
#   req(input$alternative)
# })
numK <- reactiveVal()
observe({
  req(data()$matrix)
  newValue <- ifelse(input$refFamily == "Beta", 
                     round(2*nrow(req(data()$matrix))/100),
                     nrow(req(data()$matrix)))
  numK(newValue)
})

observeEvent(input$buttonValidate, {
  newValue <- req(input$valueK)
  numK(newValue)
})
# numK <- eventReactive(input$buttonValidate, {
#   req(input$valueK)
# })
cal <- reactive({
  withProgress(value = 0, message = "Perform calibration ... ", {
  incProgress(amount = 0.3)
  t1 <- Sys.time()
  cal <- R.cache::memoizedCall(calibrateJER,
                        req(data()$matrix), categ = data()$categ, 
                        B = numB(), alpha = alpha(), 
                        refFamily = refFamily(), alternative = alternative(), 
                        K = numK()
  )
  t2 <- Sys.time()
  print(paste("calibration :",difftime(t2, t1)))
  setProgress(value = 0.7, detail = "Done")
  return(cal)
  # calibrateJER(data()$matrix, categ = data()$categ, 
  #              B = numB(), alpha = alpha(), 
  #              refFamily = refFamily(), alternative = alternative(), 
  #              K = numK()
  # )
  })
})

thr <- reactiveVal()

observe({ # if gene expression data matrix is used
  req(data()$matrix)
  newValue <- req(cal()$thr)
  thr(newValue)
})

observe({ # if p.value matrix is used
  req(data()$df)
  m = dim(data()$df)[1]
  newValue <- SimesThresholdFamily(m, kMax = m)(alpha())
  thr(newValue)
})

df <- reactiveVal()

observe({
  req(data()$matrix)
  dex <- rowWelchTests(req(data()$matrix), data()$categ)
  pval <- dex[["p.value"]]
  logp <- -log10(pval)
  fc <- dex$meanDiff
  adjp <- p.adjust(pval, method = "BH")
  newValue <- list(logp = logp, fc = fc, adjp = adjp, pval = pval)
  df(newValue)
})

observe({
  req(data()$df)
  pval <- data()$df[['p.value']]
  fc <- data()$df[['fc']]
  logp <- -log10(pval)
  adjp <- p.adjust(pval, method = "BH")
  newValue <- list(logp = logp, fc = fc, adjp = adjp, pval = pval)
  df(newValue)
})

observe({
  if(is.null(data()$df) & is.null(data()$matrix)){
    newValue <- NULL
    df(newValue)
  }
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
  FP1 <- maxFP(df()$pval[selectedGenes()$sel1], thr = thr())
  TP1 <- n1 - FP1
  FDP1 <- round(FP1/max(n1, 1), 2)
  
  n2 <- length(selectedGenes()$sel2)
  FP2 <- maxFP(df()$pval[selectedGenes()$sel2], thr = thr())
  TP2 <- n2 - FP2
  FDP2 <- round(FP2/max(n2, 1), 2)
  
  n12 <- length(selectedGenes()$sel12)
  FP12 <- maxFP(df()$pval[selectedGenes()$sel12], thr = thr())
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
    min05 <- which.min(abs(df()$adjp-i))
    y05 <- df()$logp[[min05]]
    listLog <- c(listLog, y05)
  }
  return(listLog)
})

output$lineAdjp <- renderPrint({list(lineAdjp(), class(lineAdjp()))})

thr_yaxis <- reactive({
  thrYaxis(thr = thr(), maxlogp=max(df()$logp))
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



posteriori <- 
  reactive({
    # eventReactive(input$buttonValidate, {
    
    setProgress(value = 0.9, detail = "posteriori Reactive ... ")
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
        dragmode = "select",
        showlegend = FALSE) %>%
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
  withProgress( message = "Posterio plot ...", value = 0, {
    p <- posteriori()
    shiny::setProgress(value = 1, detail = "Done")
    return(p)
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

tableCSV <- reactiveVal()
observe({
  newValue <- data.frame(row.names = base::rownames(req(data()$matrix)))
  tableCSV(newValue)
})
observeEvent(input$resetCSV,{
  tableCSV(data.frame(row.names = base::rownames(req(data()$matrix))))
})
observeEvent(d(), { 
  req(data()$matrix)
  vecteur <- rep(0, dim(data()$matrix)[1])
  vecteur[manuelSelected()] <- 1
  nameCol <- colnames(tableCSV())
  df <- cbind(tableCSV(), selection2=vecteur)
  colnames(df) <- c(nameCol, paste("User selection", length(nameCol)+1))
  tableCSV(df)
})


calcBoundSelection <- reactive({ #calculate bounds of selected genes 
  req(manuelSelected())
  calcBounds(df()$pval[manuelSelected()], thr = thr())
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
  plotMaxFP(pval = df()$pval[selectedGenes()$sel12], thr = thr()) + 
    ggtitle("Upper Left + right")
})

output$curveMaxFPSelect <- renderPlotly({
  if(length(userDTselectPost()) == 1){
    plotMaxFP(pval = df()$pval[selectionUserRe()$sel], thr = thr()) + 
      ggtitle(userDTselectPost())
  }else{
    plotMaxFP(pval = df()$pval[manuelSelected()], thr = thr()) + 
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
  req(data()$biologicalFunc)
  boundGroup(df(), 
             data()$biologicalFunc, 
             thr = thr(),
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

priori <- 
  reactive({
    # eventReactive(input$buttonValidate, {
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
  withProgress( message = "Plot", {
    priori()
    shiny::setProgress(value = 1, detail = "Done")
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
  plotMaxFP(pval = df()$pval[selectionGroup()$sel], thr = thr()) + 
    ggtitle(userDTselectPrio()) 
})



})
