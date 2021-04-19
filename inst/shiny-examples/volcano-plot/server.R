shinyServer(function(input, output, session) {
  
  source("function.R")
  
  
  # size of input data sets
  options(shiny.maxRequestSize=1024^3)
  output$help <- renderUI({
    a("IIDEA help page", href = "https://pneuvial.github.io/sanssouci/articles/IIDEA.html", target= "_blank")
  })
  
  ################### 
  # Loading data
  ###################
  
  
  ## loading example of exression gene matrix from GSEABenchmarkeR::loadEData,"geo2kegg"
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
  
  ## Loading gene set from EnrichmentBrowser::getGenesets for geo2kegg
  go.gs <- reactive({
    withProgress(message = "Load EnrichmentBrowser getGenesets ... ", {
      t1 <- Sys.time()
      go.gs <- R.cache::memoizedCall(EnrichmentBrowser::getGenesets,
                                     org = "hsa", db = "go", onto = "BP", mode = "GO.db")
      T3 <- Sys.time()
      go.gs <- R.cache::memoizedCall(cleanGo.GS, go.gs)
      t2 <- Sys.time()
      print(paste("Load EnrichmentBrowser getGenesets:",difftime(t2, t1)))
      print(paste("T3-T1:",difftime(T3, t1)))
      print(paste("T2-T3:", difftime(t2, T3)))
      setProgress(value = 1, detail = "Done")
      return(go.gs)
    })
  })
  
  ## button run 
  
  isolate({shinyjs::disable("buttonValidate")}) #while geo2kegg is not loaded, user cannot "run" #Initialization
  
  isolate(go.gs())
  
  observeEvent(geo2kegg(),{ # geo2kegg is loaded, user can "run"
    shinyjs::enable("buttonValidate")
    
  })
  
  
  ## input for example data sets
  
  nameGeo2kegg <- reactive({namedGeo2kegg(geo2kegg())}) # get names of data sets
  
  output$choiceGSEAUI <- renderUI({ #create input 
    selectInput("choiceGSEA", label = "Choose a gene data set", 
                choices = c('Leukemia (ALL): BCR/ABL mutated vs wild type'='OurData', nameGeo2kegg()))
  })
  
  ## Download our example data set, loaded from sansSouci.data
  
  ### cleaning of data set
  exampleData <- reactive({
    withProgress(value = 0, message = "Preparation for downloaded data ...", {
      
      #### expression matrix
      setProgress(value = 0.1, detail ="Download data ... ")
      data <- list()
      data$matrix <- expr_ALL
      
      setProgress(value = 0.4, detail = "Cleaning data ... ")
      
      categ <- colnames(data$matrix)
      data$categ <- rep(1, length(categ))
      data$categ[which(categ == "NEG")] <- 0
      
      colnames(data$matrix) <- data$categ
      
      data$geneNames <- base::rownames(data$matrix)
      
      setProgress(value = 0.4, detail = "Fc and p.value matrix ... ")
      
      #### degraded matrix
      dex <- rowWelchTests(data$matrix, data$categ)
      data$degrade = data.frame("p.value"=dex[["p.value"]], "fc"=dex$meanDiff)
      rm(dex)
      
      #gene set matrix
      setProgress(value = 0.7, detail = "Preparation of gene set data ...  ")
      bioFun <- expr_ALL_GO
      stopifnot(nrow(bioFun) == nrow(data$matrix))  ## sanity check: dimensions
      mm <- match(base::rownames(bioFun), base::rownames(data$matrix))
      stopifnot(!any(is.na(mm)))
      data$biologicalFunc <- bioFun[mm, ]
      rm(bioFun)
      setProgress(value = 1, detail = "Done")
      return(data)
    })
  })
  
  ### load data sets in the button : save in three csv file containing in a zip file
  output$downloadExampleData <- downloadHandler(
    filename = function() {
      paste("ExampleData", "zip", sep=".")
    },
    content = function(fname) {
      fs <- c()
      tmpdir <- tempdir()
      setwd(tempdir())
      paths <- c("expressData.csv", "biologicalFunction.csv", "degradedData.csv")
      write.csv(exampleData()$matrix, paths[1])
      write.csv(exampleData()$biologicalFunc, paths[2])
      write.csv(exampleData()$degrade, paths[3])
      zip(zipfile=fname, files=paths)
    },
    contentType = "application/zip"
  )
  
  ## Loading user inputs
  
  ### Loading expression gene matrix
  fileData <- reactiveVal(NULL) #Initialisation
  observe({ #when a new csv file is in input
    newValue <- req(input$fileData)
    fileData(newValue)
  })
  observeEvent(input$resetInputData, { #to delete input file (clicking bin incon)
    fileData(NULL) # delete serveur variable
    reset("fileData") #delete UI variable (in fileInput)
  })
  
  ## Loading gene set matrix : same structure than before
  fileGroup <- reactiveVal(NULL)
  observe({
    newValue <- input$fileGroup
    fileGroup(newValue)
  })
  observeEvent(input$resetInputGroup, {
    fileGroup(NULL)
    reset("fileGroup")
  })
  
  ## Cleaning data
  data <- 
    eventReactive(input$buttonValidate, {
      withProgress(value = 0, message = "Upload Data... ", {
        data <- list()
        if (input$checkboxDemo){ #if example data set
          if(req(input$choiceGSEA)=='OurData'){ # cleaning for data from sanssouci.data
            setProgress(value = 0.4, detail = "SansSoucis data set ...")
            data$matrix <- expr_ALL #read data from sansSouci.data
            
            ### cleaning categories 
            categ <- colnames(data$matrix)
            data$categ <- rep(1, length(categ))
            data$categ[which(categ == "NEG")] <- 0
            
            data$geneNames <- rownames(data$matrix)
            
            setProgress(value = 0.7, detail = "Preparation of gene set data ...  ")
            
            bioFun <- expr_ALL_GO
            stopifnot(nrow(bioFun) == nrow(data$matrix))  ## sanity check: dimensions
            mm <- match(base::rownames(bioFun), base::rownames(data$matrix))
            stopifnot(!any(is.na(mm)))
            data$biologicalFunc <- bioFun[mm, ]
            rm(bioFun)
            data$boolValidation <- TRUE 
          } else { # cleaning data set from GSEA data set
            setProgress(value = 0.4, detail = "GSEA data set ...")
            rawData <- R.cache::memoizedCall(maPreproc,geo2kegg()[input$choiceGSEA])[[1]]
            
            data$matrix <- SummarizedExperiment::assays(rawData)$exprs
            
            cats <- SummarizedExperiment::colData(rawData)
            ww <- match(cats$Sample, base::colnames(data$matrix))
            categ <- cats$GROUP[ww]
            data$categ <- categ
            setProgress(value = 0.7, detail = "GSEA data set ...")
            
            data$geneNames <- base::rownames(data$matrix)
            
            data$biologicalFunc <- go.gs() #On a laissé sous forme de liste car on a adapté les fonctions qui en ont besoin. Plus rapide qu'en la transformant en matrice binaire
            
            data$boolValidation <- TRUE
          }
          
        } else { # if user use his own data
          req(input$fileData)
          req(fileData())
          file <- req(fileData())
          setProgress(value = 0.1, detail = "Read csv ...")
          data$matrix <- read.csv(file = file$datapath, row.names = 1, check.names=FALSE)
          
          data$boolDegrade <- (length(data$matrix) == 2 & all(sort(colnames(data$matrix)) == sort(c("fc","p.value")))) #est ce que la matrice est dégradée
          setProgress(value = 0.2, detail = "Test data ...")
          if (data$boolDegrade){ #cleaning degraded matrix
            setProgress(value = 0.4)
            data$df <- data$matrix
            data$geneNames <- base::rownames(data$matrix)
          } else { #cleaning expression matrix
            setProgress(value = 0.4, detail = "Clean data set ...")
            clean <- cleanMatrix(data$matrix) #see cleanMatrix function : 
            ### gve specific issue for non available matrix
            
            ### if matrix is not available, data$matrix == NULL allows to block the calibration JER and the printing  
            if(clean$boolValidation){ #if matrix is ok
              data$matrix = clean$data
            }else{ #if isn't 
              data$matrix=NULL
            }
            
            data$matrix.color <- clean$color #color of error message
            data$matrix.text <- clean$text #eroor message
            
            data$categ <- colnames(data$matrix)
            
            
            data$geneNames <- base::rownames(data$matrix)
            data$boolValidation <- clean$boolValidation
          }
          
          setProgress(value = 0.7, detail = "Preparation of gene set data ...")
          
          ## cleaning of gene set matrix
          fileGroup <- fileGroup()
          if (!is.null(fileGroup)){
            setProgress(value = 0.75, detail = "Read gene set data ...")
            T1 <- Sys.time()
            bioFun <- read.csv(file = fileGroup$datapath, row.names = 1, check.names=FALSE)
            T2 <- Sys.time()
            print(paste("read gene set", T2-T1))
            setProgress(value = 0.8, detail = "Cleaning of gene set data ...")
            cleanBio <- cleanBiofun(bioFun) # cleaning and message error
            
            data$bioFun.color <- cleanBio$color
            data$bioFun.text <- cleanBio$text
            
            bioFun <- cleanBio$biofun
            setProgress(value = 0.9, detail = "Matching ...")
            
            matchBio <- matchMatrixBiofun(matrixFunc = data$matrix, biofun = bioFun) #verification compatibility between the two matrices
            if(matchBio$boolValidation & cleanBio$boolValidation){ #if ok
              data$biologicalFunc <- matchBio$biofun
            } else { #if not ok
              data$biologicalFunc <- NULL
            }
            
            data$match.color <- matchBio$color
            data$match.text <- matchBio$text
            
            rm(bioFun)
          }
          
          # if matrix is not available (here degraded matrix), data$matrix == NULL allows to block the calibration JER and the thr is calculate with simes and k=m 
          if (data$boolDegrade){
            data$matrix <- NULL 
          } 
          
          
          
        }
        setProgress(value = 1, detail = "Done")
        return(data)
      })
    }
    )
  
  matrixChosen <- reactive ({ ## to keep nrow of the chosen matrix
    boolDegrade <- FALSE
    if (input$checkboxDemo){
      if(req(input$choiceGSEA)=='OurData'){
        matrix <- expr_ALL
        
      } else {
        rawData <- R.cache::memoizedCall(maPreproc,geo2kegg()[input$choiceGSEA])[[1]]
        matrix <- SummarizedExperiment::assays(rawData)$exprs
      }
      
    } else { # if user use his own data
      req(input$fileData)
      req(fileData())
      file <- req(fileData())
      matrix <- read.csv(file = file$datapath, row.names = 1, check.names=FALSE)
      boolDegrade <- (length(matrix) == 2 & all(sort(colnames(matrix)) == sort(c("fc","p.value")))) #est ce que la matrice est dégradée
    }
    return(list(matrix = matrix, boolDegrade = boolDegrade))
  })
  
  
  urlDataSet <- eventReactive(input$buttonValidate, { #give link to description of geo2kegg data sets
    req(geo2kegg())
    req(input$choiceGSEA)
    req(geo2kegg()[[input$choiceGSEA]])
    return(a("URL link to data set description", href=geo2kegg()[[input$choiceGSEA]]@experimentData@url))
    
  })
  output$msgURLds <- renderUI({ 
    req(urlDataSet)
    tagList(urlDataSet())
  })
  
  ## different error message for non compliant matrix
  output$errorInput <- renderUI({ # express gene matrix error
    tags$span(style= req(data()$matrix.color), paste(req(data()$matrix.text)))
  })
  
  output$errorBioMatrix <- renderUI({ #gene set matrix error
    tags$span(style=req(data()$bioFun.color), paste(req(data()$bioFun.text)))
  })
  output$errorMatch <- renderUI({ # none commun gene btw expression gene matrix and gene set matrix
    tags$span(style=req(data()$match.color), paste(req(data()$match.text)))
  })
  
  
  
  
  ################### 
  # Parameters 
  ###################
  
  # Confidance alpha
  alpha <- reactiveVal(0.1) #Initialization
  observeEvent(input$buttonValidate,{ # When Run is clicked : we get the input value
    newValue <- req(1 - input$sliderConfLevel/100)
    alpha(newValue)
  })
  
  # number of permutation
  numB <- reactiveVal(500) #Initialization
  observeEvent(input$buttonValidate, { # When Run is clicked : we get the input value
    newValue <- req(input$numB)
    numB(newValue)
  })
  
  # reference family 
  refFamily <- reactiveVal("Simes") #Initialization
  observeEvent(input$buttonValidate, { # When Run is clicked : we get the input value
    newValue <- req(input$refFamily)
    refFamily(newValue)
  })
  
  #alternative hypothesis
  alternative <- reactiveVal("two.sided") #Initialisation 
  observeEvent(input$buttonValidate, { # When Run is clicked : we get the input value
    newValue <- req(input$alternative)
    alternative(newValue)
  })
  
  # parameter K
  ## dynamic input (need matrixChosen())
  output$inputK <- renderUI({
    req(matrixChosen())
    numericInput("valueK", 
                 label = "K (size of reference family)", 
                 value = numKI(),
                 min = 1,
                 max = nrow(matrixChosen()$matrix))
  })
  
  ## numKI() is used to intiate the printed input valueK 
  numKI <- reactiveVal()
  observe({ #Initialization, if refFamily == 'Beta' or not 
    req(matrixChosen())
    newValue <- ifelse(input$refFamily == "Beta", 
                       round(2*req(nrow(matrixChosen()$matrix))/100),
                       req(nrow(matrixChosen()$matrix)))
    numKI(newValue)
  })
  
  ## numK is the parameters choosen by users and use in server side
  numK <- reactiveVal()
  observeEvent(data(),{   #si les paramètres ne sont pas activés, input$valueK n'existe pas et donc on a pas le bon résultat ... pourquoi ? CalibrateJER devrait gérer si K est null non ? 
    req(data())           # when data() is change, INITIALISATION of numK()
    newValue <- req(nrow(data()$matrix))
    numK(newValue)
  })
  observeEvent(input$buttonValidate, { # When Run is clicked : we get the input value. Here an issue : if, advanced parameters is not opened, input$valueK == NULL
    newValue <- req(input$valueK)
    numK(newValue)
  })
  
  
  # If degraded matrix is available 
  
  output$msgDegraded <- renderUI({
    tags$span(style= "color:grey", paste("A matrix containing p-values and fold change is detected.", 
                                         "Thus, you cannot change the following advanced parameters:\n",
                                         "Reference family = 'Simes',\n K = ", nrow(matrixChosen()$matrix)))
  })
  
  observe({
    req(matrixChosen())
    if(matrixChosen()$boolDegrade){
      shinyjs::show("msgDegraded")
      
      shinyjs::hide("alternative")
      shinyjs::hide("numB")
      shinyjs::hide("refFamily")
      shinyjs::hide("inputK")
    } else {
      shinyjs::hide("msgDegraded")
      
      shinyjs::show("alternative")
      shinyjs::show("numB")
      shinyjs::show("refFamily")
      shinyjs::show("inputK")
    }
  })
  
  
  ################### 
  # Calibration
  ###################
  
  
  
  ## JER calibration on available expression gene matrix
  cal <- reactive({
    # eventReactive(input$buttonValidate, {
    withProgress(value = 0, message = "Perform calibration ... ", {
      incProgress(amount = 0.3)
      t1 <- Sys.time()
      # print(paste("alpha ", alpha()))
      # print(paste('b = ', numB()))
      # print(paste("alternative = ", alternative()))
      # print(paste("ref framily = ", refFamily()))
      # print(paste("cal : k = ", numK()))
      cal <- R.cache::memoizedCall(calibrateJER,
                                   req(data()$matrix), # if data()$matrix == NULL, not perform [-> not available matrix or degraded matrix]
                                   categ = data()$categ, 
                                   B = numB(), alpha = alpha(), 
                                   refFamily = refFamily(), alternative = alternative(), 
                                   K = numK()
      )
      t2 <- Sys.time()
      print(paste("calibration :",difftime(t2, t1)))
      setProgress(value = 0.7, detail = "Done")
      return(cal)
    })
  })
  
  ## Thr calculation
  thr <- reactiveVal() #Initialization
  
  observe({ 
    req(data()$matrix) # if gene expression matrix is available 
    newValue <- req(cal()$thr) #we use thr from calibration and user's parameters
    thr(newValue)
  })
  
  observe({ # if p.value matrix is used
    req(data()$df) #if degraded matrix is used
    # if the gene expression matrix is available => data()$df == NULL 
    m = dim(data()$df)[1]
    newValue <- SimesThresholdFamily(m, kMax = m)(alpha()) # force using of Simes and k=m
    thr(newValue)
  })
  
  ################### 
  # P-values and foldchange
  ###################
  
  df <- reactiveVal() #creation
  
  observe({ # If gene expression matrix is available
    req(data()$matrix) #active observe cell
    dex <- rowWelchTests(req(data()$matrix), data()$categ) # make a test
    pval <- dex[["p.value"]]
    logp <- -log10(pval)
    fc <- dex$meanDiff
    adjp <- p.adjust(pval, method = "BH")
    newValue <- list(logp = logp, fc = fc, adjp = adjp, pval = pval)
    df(newValue) #update df()
  })
  
  observe({ # if degraded version
    req(data()$df) #active observe cell
    pval <- data()$df[['p.value']] # from input
    fc <- data()$df[['fc']]
    logp <- -log10(pval)
    adjp <- p.adjust(pval, method = "BH")
    newValue <- list(logp = logp, fc = fc, adjp = adjp, pval = pval)
    df(newValue)
  })
  
  observe({ #if none matrix is in input
    if(is.null(data()$df) & is.null(data()$matrix)){
      newValue <- NULL
      df(newValue)
    }
  })
  
  ################### 
  # Threshold 
  ###################
  
  # vertical contient la valeur de déplacement de tous les objets déplaçable (dans notre cas juste les seuils)
  vertical <- reactive({event_data("plotly_relayout", source='A')})
  
  #####
  # Contenu de vertical : 
  # vertical()[["shapes[0].x0"]] : seuil logFC droit
  # vertical()[["shapes[2].x0"]] : seuil logFC gauche 
  # vertical()[["shapes[1].y0"]] : seuil pvalue 
  #
  # vertical()[["shapes[A].y0"]] : A pour le numéro de l'objet et x0 ou y0 pour le sens de déplacement
  
  ## threshold logfc right
  xint <- reactiveVal(0.5)
  observeEvent(vertical()[["shapes[0].x0"]], { #activate when the object 0 change
    if(input$symetric){ 
      newValue <- vertical()[["shapes[0].x0"]]
      xint(newValue)
      xint2(-newValue)
    } else {
      newValue <- vertical()[["shapes[0].x0"]]   
      xint(newValue)  
    }
  })
  
  ## threshold logfc left
  xint2 <- reactiveVal(-0.5)
  observeEvent(vertical()[["shapes[2].x0"]], {#activate when the object 2 change
    if( input$symetric){
      newValue <-  vertical()[["shapes[2].x0"]]
      xint( - newValue)
      xint2(newValue)
    } else {
      newValue <- vertical()[["shapes[2].x0"]] 
      xint2(newValue)
    }
    
  })
  
  ## threshold pval
  yint <- reactiveVal()
  observeEvent(df(), { # initialisation : adjusted p-value == 0.1
    min0.1 <- which.min(abs(df()$adjp-0.1))
    y0.1 <- df()$logp[[min0.1]]
    print(paste("min0.1 = ",min0.1," y0.1 = ", y0.1))
    yint(y0.1)
  })
  observeEvent(vertical()[["shapes[1].y0"]], { # when user change threshold on plotly
    p <- 10^(-vertical()[["shapes[1].y0"]])
    y_sel <- which((df()$pval <= p))          ## selected by  p-value
    newValue <- Inf
    if (length(y_sel) > 0) {
      newValue <- min(df()$logp[y_sel])              ## threshold on the log(p-value) scale
    }
    yint(newValue)
  })
  
  ################### 
  # Selecting gene with thresholds
  ###################
  
  # selected genes for server calcuation
  selectedGenes <- reactive({
    req(df())
    ## gene selections
    sel1 <- which(df()$logp >= yint() & df()$fc >= xint()) #upper right 
    sel2 <- which(df()$logp >= yint() & df()$fc <= xint2()) #upper left 
    sel12 <- sort(union(sel1,sel2)) #both
    return(list(sel1 = sel1, sel2 = sel2, sel12 = sel12))
  })
  
  
  #matrix p-val & fc of selected genes by thresholds => to plot red points on volcano plot
  selected_points <- reactive({
    list(x = df()$fc[selectedGenes()$sel12], 
         y = df()$logp[selectedGenes()$sel12])
  })
  
  ################### 
  # Selecting gene by lasso or box 
  ###################
  
  # reactive variable containing value of lasso/box selecting
  d <- reactive({ event_data("plotly_selected", source='A')})
  
  # list of selected genes (list of name rows)
  manuelSelected <- reactive({
    req(d())
    d()[['pointNumber']][which(d()[['curveNumber']]==0)]+1
  })
  
  
  ################### 
  # Calculate post hob bounds
  ###################
  
  # POST Hoc bound (PHB) for thresholds selection
  TP_FDP <- reactive({
    req(selectedGenes())
    ## post hoc bounds in selections 
    ### upper right
    # n1 <- length(selectedGenes()$sel1) # number of selected genes
    # FP1 <- maxFP(df()$pval[selectedGenes()$sel1], thr = thr()) # 
    # TP1 <- n1 - FP1 #
    # FDP1 <- round(FP1/max(n1, 1), 2) #
    ### upper left
    # n2 <- length(selectedGenes()$sel2)
    # FP2 <- maxFP(df()$pval[selectedGenes()$sel2], thr = thr())
    # TP2 <- n2 - FP2
    # FDP2 <- round(FP2/max(n2, 1), 2)
    ### both
    n12 <- length(selectedGenes()$sel12)
    FP12 <- maxFP(df()$pval[selectedGenes()$sel12], thr = thr())
    TP12 <- n12 - FP12
    FDP12 <- round(FP12/max(n12, 1), 2)
    
    return(list(
      # n1 = n1, TP1 = TP1, FDP1 = FDP1,
      # n2 = n2, TP2 = TP2, FDP2 = FDP2,
      n12 = n12, TP12 = TP12, FDP12 = FDP12))
    
  })
  
  # calculate PHB for lasso box selection
  calcBoundSelection <- reactive({ # 
    req(manuelSelected())
    calcBounds(df()$pval[manuelSelected()], thr = thr()) #return(list(n=n, FP=FP, TP=TP, FDP=FDP))
  })
  
  
  
  ################### 
  # Post Hoc bound table reactive variable
  ###################
  
  #Initialization of variable of post hoc bound table (PHB table)
  tableResult <- reactiveVal(data.frame(
    Selection = c("Threshold selection"))) #Initialization
  
  # PHB table for only selected genes by thresholds
  ## reactive : TP_FDP() <- selectedGenes() <- {xint(), xint2(), yint()} <- vertical() <- user threshold moving
  baseTable <- reactive({
    data.frame(`Selection` = c("Threshold selection"), 
               "# genes" = c(TP_FDP()$n12),
               "TP≥" = as.integer(c(TP_FDP()$TP12)), 
               "FDP≤" = c(TP_FDP()$FDP12),
               check.names = FALSE)
  })
  
  # updating PHB table when thresholds change
  observeEvent(TP_FDP(),{ # When threshold change
    
    bottomTable <- tableResult() %>%  #keep box/lasso selection
      # filter(Selection != "Upper right") %>% 
      # filter(Selection != "Upper left") %>% 
      filter(Selection != "Threshold selection") #remove row named "Threshold selection"
    upperTable <- baseTable() # take new value of PHB for the new threshold selection
    newValue <- rbind(upperTable, bottomTable)
    tableResult(newValue)
  })
  
  # Add PHB for lasso and box selection 
  ## d() is the reactive variable for box and lasso selection
  observeEvent(d(),{  # When user selects a new group of points
    req(calcBoundSelection())
    vectorGene <- names(df()$pval[manuelSelected()]) #list of gene contained in gene selection
    url <- UrlStringdbGrah(vectorGene) #construction of link to StringDB interaction graph
    n <- dim(tableResult())[1]
    newValue <- rbind(tableResult(), c(paste('<a target="_blank" href="', url, '" >User selection ',n, '</a>', sep=""), 
                                       calcBoundSelection()$n, 
                                       calcBoundSelection()$TP,
                                       round(calcBoundSelection()$FDP, 2)))
    tableResult(newValue)
  })
  
  # To clean gene selection from lasso/box selection 
  ## keep only threshold selection
  observeEvent(input$resetCSV, { # to clean printed table
    newValue <- baseTable()
    tableResult(newValue)
  })
  
  # if data change, PHB table is cleaned
  observeEvent(data(), { # to clean printed table
    newValue <- baseTable()
    tableResult(newValue)
  })
  

  ################### 
  # Post Hoc bound table outputs
  ###################
  
  # reactive popify to explain PHB table
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
  
  # output for PHB table 
  output$tableBounds <- renderDT({
    req(TP_FDP())
    
    tableResult()
    
  }, selection = list(mode = 'single', selectable = -(1))  # can select only one row and not the first one
  , escape = FALSE # to print url link to open stringDB graph
  )
  
  ################### 
  # Prepare objects for plotly volcano plots
  ###################

  #value for adjusted p-values yaxis
  lineAdjp <- reactive({ 
    req(df())
    listLog <- c()
    for (i in c(0.5,0.25,0.1,0.05,0.025,0.01,0.001,0.0001)){ #selected line on yaxis /!\ if you change it you have to change in yaxis()
      min05 <- which.min(abs(df()$adjp-i)) #to take the gene with the nearest adjpvalue to i 
      y05 <- df()$logp[[min05]] # take the equivalent in logpvalue to print it on VP
      listLog <- c(listLog, y05)
    }
    return(listLog)
  })
  
  # value for 'NUmber of false discoverie' yaxis : optimize y line
  thr_yaxis <- reactive({
    req(alpha()) # if parameters change (not only alpha)
    req(thr())
    req(df())
    thrYaxis(thr = thr(), maxlogp=max(df()$logp))
  })
  
  #reactive variable containing values for yaxis depending users choice (input$choiceYaxis)
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
                      ticktext = c(0.5,0.25,0.1,0.05,0.025,0.01,0.001,0.0001) #be careful of changing of lingAdjp()
                    ),
                    "thr" = list(
                      title = "Maximal number of false positives", 
                      titlefont = f, 
                      autotick = FALSE,
                      tickmode = "array",
                      tickvals = thr_yaxis()$pvalue,
                      ticktext = thr_yaxis()$num
                    )
    )
    return(yaxis)
  })
  
  # reactive values for threshold (used for selecting genes)
  thrLine <- reactive({
    list(
      list( # rigth logFC threshold vertical()[["shapes[0].x0"]]
        type = "line",
        line = list(color = "orange", dash = "dot"),
        x0 = xint(),
        x1 = xint(),
        y0 = 0,
        y1 = 1,
        yref = "paper"
      )
      ,
      list( # p-val threshold vertical()[["shapes[1].y0"]]
        type = "line",
        line = list(color = "orange", dash = "dot"),
        x0 = 0,
        x1 = 1,
        y0 = yint(),
        y1 = yint(),
        xref = "paper"
      )
      ,
      list( # left logFC threshold vertical()[["shapes[2].x0"]]
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
  
  ################### 
  # OUTPUT object of plotly volcano plot VP1
  ###################
  
  # Intialisation of volcano plot for threshold selection (part 1)
  posteriori <- reactive({
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
                            color = 'grey'), 
              name = 'unselected',
              type='scattergl', mode = "markers", #make points
              source='A' #name plotly::plot
              # , text = annotation()[['nameGene']],
              # customdata = paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", annotation()[['nameGene']])
              # , height = 600
      ) %>%
        add_markers(x = selected_points()$x, y = selected_points()$y, #add at the begening selected points in red
                    marker = list(
                      color = "red",
                      size = 6
                    ),
                    name ="selected") %>%
        layout(
          xaxis = list(title = "Fold change (log scale)", titlefont = f),
          yaxis = isolate(yaxis()), #isolate is used not to reactive this reactive variable (improve global reactivity : see bellow for the changing of yaxis)
          title = "",
          shapes = isolate(thrLine()), #same as before
          dragmode = "select", #initialise lasso/box selection by default
          showlegend = FALSE) %>%
        # onRender("
        #           function(el) {
        #               el.on('plotly_click', function(d) {
        #                   var url = d.points[0].customdata;
        #                   window.open(url);
        #               });
        #           }
        # ") %>% #active clicking on points
        event_register("plotly_selecting") %>% # to save point selected by lasso:box selection
        config(editable = TRUE)%>%
        toWebGL() #to go faster on web
      
      
    })
  
  

  output$volcanoplotPosteriori <- renderPlotly({
    withProgress( message = "Posterio plot ...", value = 0, {
      p <- posteriori()
      shiny::setProgress(value = 1, detail = "Done")
      return(p)
    })
  })
  
  
  
  ################### 
  # Change of VP1
  ###################
  
  # This part aims of changing only hight layers of the VP1 plotly object. 
  # Indeed, posteriori() is only reactive if df() change. 
  # Each 'observeEvent' represents an user change. 
  # The VP1 is updated with the function plorlyProxy (choose the plotly OUTPUT) 
  # and the plotlyProxyInvoke() function to define the action.
  
  
  #when threshold moved and red selected points change
  observeEvent( vertical(), { 
    
    plotlyProxy("volcanoplotPosteriori", session) %>%
      plotlyProxyInvoke("deleteTraces",1)%>% #delete the previous red points
      plotlyProxyInvoke(
        "addTraces", #add red points from selected_points()
        list(
          x = unname(selected_points()$x),
          y = unname(selected_points()$y),
          type = "scatter",
          mode = "markers",
          line = list(color = "red"), 
          name = "selected", 
          showlegend = TRUE
        ), 1 #the rank on the stack
      )
  })
  
  # Yaxis change 
  observeEvent(
    {input$choiceYaxis #user selection
      yaxis()} #update yaxis values when parameters change (alpha, ...)
    , { #when we choose a different y axis 
      plotlyProxy("volcanoplotPosteriori", session) %>%
        plotlyProxyInvoke("relayout", list(yaxis = yaxis())) #relayout to update values
      #a stack is not used in layout, only one value is possible to yaxis
    })
  
  #When user changes orange threshold to select red points
  observeEvent(vertical(), {  # for moving of orange threshold
    plotlyProxy("volcanoplotPosteriori", session)%>%
      plotlyProxyInvoke("relayout", 
                        list(shapes = thrLine())) #same as before
  })
  
  # on the stack TRACES of VP1, the layer 0 is all the point (unselected), 
  #the layer 1 is the red points
  # the layer 2 is the blue points
  
  
  
  # when user want to watch again a user selection by box or lasso
  #tableBounds_rows_selected : the user selection on the output tableBounds
  ## take the name of the user selection 
  userDTselectPost <- reactive({
    href <- tableResult()[input$tableBounds_rows_selected, "Selection"]
    str_remove_all(str_remove_all(href, "<a(.*?)>"), "(</a>)")
  })
  
  #take the list of gene selected by user 
  selectionUserRe <- reactive({
    vect <- tableCSV()[,userDTselectPost()] 
    sel <- which(vect == 1)
    list(sel = sel)
  })
  
  # update plotly VP1
  observeEvent(userDTselectPost(), {
    if(length(userDTselectPost()) == 1){ #if a row is selected (length == 1 because can select only one row)
      plotlyProxy("volcanoplotPosteriori", session) %>% 
        plotlyProxyInvoke("deleteTraces", 2)            #first we delete previous blue points : 
      # if there are not bleu points, none points are deleted
      
      plotlyProxy("volcanoplotPosteriori", session) %>% #second, we trace new blue points
        plotlyProxyInvoke(
          "addTraces",
          list(
            x = unname(df()$fc[selectionUserRe()$sel]),
            y = unname(df()$logp[selectionUserRe()$sel]),
            type = "scattergl",
            mode = "markers",
            line = list(color = "blue"),
            name = userDTselectPost()
          ), 2 #put it in layer 2
        )
    } else { # If none, two or more row are selected
      plotlyProxy("volcanoplotPosteriori", session) %>%
        plotlyProxyInvoke("deleteTraces",2)
    }
  })
  
  
  ################### 
  # Download list of selected gene from box/lasso
  ###################
  
  tableCSV <- reactiveVal() #creation of reactive variable
  
  #initialize table when data() change
  observe({
    newValue <- data.frame(row.names = base::rownames(req(data()$matrix))) #a dataframe with rowname without features
    tableCSV(newValue)
  })
  
  #When user want to 
  observeEvent(input$resetCSV,{
    tableCSV(data.frame(row.names = base::rownames(req(data()$matrix))))
  })
  
  #When user selected ne gene set from lasso/box [d() activate]
  observeEvent(d(), { 
    req(data()$matrix)
    vecteur <- rep(0, dim(data()$matrix)[1])
    vecteur[manuelSelected()] <- 1
    nameCol <- colnames(tableCSV())
    df <- cbind(tableCSV(), selection2=vecteur)
    colnames(df) <- c(nameCol, paste("User selection", length(nameCol)+1))
    tableCSV(df)
  })
  
  
  
  # download binary matrix containing list of gene in user selections
  output$downloadData <- downloadHandler( #download csv of user selection
    filename = function() {
      tag <- format(Sys.time(), "%Y-%M-%d_%H-%m-%S")
      sprintf("volcano-plot_gene-selections_%s.csv", tag)
    },
    content = function(file) {
      write.csv(tableCSV(), file)
    }
  )
  
  # download PHB table
  output$downloadPHBTable <- downloadHandler( 
    filename = function() {
      tag <- format(Sys.time(), "%Y-%M-%d_%H-%m-%S")
      sprintf("volcano-plot_bounds_%s.csv", tag)
    },
    content = function(file) {
      write.csv(tableResult(), file)
    }
  )
  
  # while buttonValidate is not clicke, these button are hidden
  observeEvent(input$buttonValidate, {
    shinyjs::show("downloadPHBTable")
    shinyjs::show("resetCSV")
    shinyjs::show("downloadData")
  })
  
  # output$curveMaxFPBoth <- renderPlotly({
  #   plotMaxFP(pval = df()$pval[selectedGenes()$sel12], thr = thr()) + 
  #     ggtitle("Upper Left + right")
  # })
  # 
  # output$curveMaxFPSelect <- renderPlotly({
  #   if(length(userDTselectPost()) == 1){
  #     plotMaxFP(pval = df()$pval[selectionUserRe()$sel], thr = thr()) + 
  #       ggtitle(userDTselectPost())
  #   }else{
  #     plotMaxFP(pval = df()$pval[manuelSelected()], thr = thr()) + 
  #       ggtitle("User selection")
  #   }
  # })
  
  
  
  ################### 
  # PART 2 : gene set analyses
  ###################
  
  # PHB for all gene sets
  tableBoundsGroup <- reactive({
    withProgress(message = "tableBoundsGroup", {
      T1 <- Sys.time()
      req(data())
      req(data()$biologicalFunc)
      table <- boundGroup(df(), 
                          data()$biologicalFunc, 
                          thr = thr()
      )
      T2 <- Sys.time()
      print(paste("post hoc bounds on gene set:",difftime(T2, T1)))
    })
    return(table)
  })
  
  # calculate vounds for all features to compare with each gene sets (competitive methods)
  boundsW <- reactive({ #calcul bounds for all features
    req(df())
    calcBounds(listPval = df()$pval, thr = thr())
  })
  
  
  # dowload csv file containing PHB table of gene sets
  output$downloadPHBTableGroup <- downloadHandler( #download csv of user selection
    filename = function() {
      tag <- format(Sys.time(), "%Y-%M-%d_%H-%m-%S")
      sprintf("gene-set_bounds_%s.csv", tag)
    },
    content = function(file) {
      write.csv(tableBoundsGroup(), file)
    }
  )
  
  observe({ # show the download button when 
    req(tableBoundsGroup())
    shinyjs::show("downloadPHBTableGroup")
  })
  
  # If user choose SEA alternative
  filteredTableBoundsGroup  <- reactive({
    req(input$buttonSEA)
    if (input$buttonSEA == "competitive"){
      table <- tableBoundsGroup()
      sel <- which(table[["FDP≤"]] < boundsW()$FDP)
      newValue <- table[sel,]
      return(newValue)
    } else if (input$buttonSEA == "self"){
      table <- tableBoundsGroup()
      sel <- which(table[["TP≥"]] > 0)
      return(table[sel,])
    } else {
      return(tableBoundsGroup())
    }
  })
  
  # reactive popify to explain PHB table
  output$OutQtableBoundsGroup <- renderUI({
    req(tableBoundsGroup())
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
  
  # formatRound(DT::datatable(genesInputGvsG(),options=list(pageLength=10)),columns=c(2,3,4,5,6,7),digits=3)
  output$tableBoundsGroup <- renderDT({
    table <- filteredTableBoundsGroup()
    table[["FDP≤"]] <- round(table[["FDP≤"]], 2)
    table
    
  }, selection = 'single' , escape = FALSE , options = list(scrollX = TRUE))
  
  # name of gene set selected by user  
  userDTselectPrio <- reactive({
    req(filteredTableBoundsGroup())
    req(input$tableBoundsGroup_rows_selected)
    href <- filteredTableBoundsGroup()[input$tableBoundsGroup_rows_selected,"Name"]
    name <- str_remove_all(str_remove_all(href, "<a(.*?)>"), "(</a>)")
    return(name)
  })
  
  #list of genes selected by user
  selectionGroup <- reactive({
    req(data())
    req(df())
    
    group <- req(userDTselectPrio())
    bioFun <- data()$biologicalFunc
    if (class(bioFun)=="list"){
      ids <- bioFun[[group]]
    }else{
      ids <- which(bioFun[, group] == 1)
    }
    list(sel = ids)
  })
  
  # 'reactive" plot : as the plot before, this one should be activate once. See posteriori() for details
  #VP2
  priori <- reactive({
      f <- list(
        size = 14,
        color = "#000000"
      )
      lte <- "≤"
      gte <- "≥"
      plot_ly(data.frame(x = df()$fc, y=df()$logp), 
              x = ~x, y = ~y, 
              marker = list(size = 2,
                            color = 'grey'), 
              name = 'genes',
              type='scattergl', mode = "markers", source='B'
              ,
              text = data()$geneNames,
              customdata = paste0("http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=", data()$geneNames)
      )%>% 
        layout(
          showlegend = TRUE,
          xaxis = list(title = "Fold change (log scale)", titlefont = f),
          yaxis = isolate(yaxis()),
          title = "",
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
  
  #output of priori()
  output$volcanoplotPriori <- renderPlotly({
    withProgress( message = "Plot", {
      p <- priori()
      shiny::setProgress(value = 1, detail = "Done")
      return(p)
    })
  })
  
  #when yaxis change
  observeEvent({input$choiceYaxis
    yaxis()}, { #when we choose a different y axis 
      plotlyProxy("volcanoplotPriori", session) %>%
        plotlyProxyInvoke("relayout", list(yaxis = yaxis()))
      
      
    })
  
  
  #when user select a gen set to print it on VP2
  #here, there are not red points : stack is composed of 0 :points ; 1 : blue points (gene set selected)
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
  
  # output$curveMaxFPGroup <- renderPlotly({
  #   req(selectionGroup())
  #   plotMaxFP(pval = df()$pval[selectionGroup()$sel], thr = thr()) + 
  #     ggtitle(userDTselectPrio()) 
  # })
  
  
  
})
