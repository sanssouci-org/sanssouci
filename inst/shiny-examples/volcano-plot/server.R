shinyServer(function(input, output, session) {
  
  source("function.R")
  
  
  # size of input data sets
  options(shiny.maxRequestSize=1024^3)
  output$help <- renderUI({
    a("IIDEA help page", href = "https://sanssouci-org.github.io/sanssouci/articles/IIDEA.html", target= "_blank")
  })
  
  ################### 
  # Loading data
  ###################
  
  
  ## loading example of exression gene matrix from GSEABenchmarkeR::loadEData,"geo2kegg"
  # geo2kegg <- reactive({
  #   withProgress(message = "Load GSEABenchmarkeR data set ... ", {
  #     t1 <- Sys.time()
  #     data <- R.cache::memoizedCall(GSEABenchmarkeR::loadEData,"geo2kegg")
  #     
  #     t2 <- Sys.time()
  #     print(paste("Load GSEABenchmarkeR data set:",difftime(t2, t1)))
  #     setProgress(value = 1, detail = "Done")
  #     return(data)
  #   })
  # })
  
  ## Loading gene set from EnrichmentBrowser::getGenesets for geo2kegg
  # go.gs <- reactive({
  #   withProgress(message = "Load EnrichmentBrowser getGenesets ... ", {
  #     t1 <- Sys.time()
  #     go.gs <- R.cache::memoizedCall(EnrichmentBrowser::getGenesets,
  #                                    org = "hsa", db = "go", onto = "BP", mode = "GO.db")
  #     T3 <- Sys.time()
  #     go.gs <- R.cache::memoizedCall(cleanGo.GS, go.gs) # our func
  #     t2 <- Sys.time()
  #     print(paste("Load EnrichmentBrowser getGenesets:",difftime(t2, t1)))
  #     print(paste("T3-T1:",difftime(T3, t1)))
  #     print(paste("T2-T3:", difftime(t2, T3)))
  #     setProgress(value = 1, detail = "Done")
  #     return(go.gs)
  #   })
  # })
  
  ## button run 
  
  isolate({shinyjs::disable("buttonValidate")}) #while geo2kegg is not loaded, user cannot "run" #Initialization
  
  # isolate(go.gs())
  
  # observeEvent(geo2kegg(),{ # geo2kegg is loaded, user can "run"
  observeEvent(object_I(),{ # object is loaded, user can "run"
    shinyjs::enable("buttonValidate")
    
  })
  
  
  ## input for example data sets
  
  namesExampleFile <- reactive({
    filenames <- (list.files("GSEABenchmarkeR/express-data-set", pattern="*.RDS", full.names=TRUE))
    if(length(filenames) == 0){
      return(NULL)
    }
    ldf <- lapply(filenames[2:length(filenames)], readRDS)
    lID <- sapply(ldf,function(l){l@metadata$dataId})
    names(lID) <- paste(sapply(ldf,function(l){l@metadata$experimentData@other$disease})," (", (lID),")", sep="")
    return(lID)
  }) # get names of data sets
  
  output$choiceGSEAUI <- renderUI({ #create input 
    selectInput("choiceGSEA", label = "Choose a gene data set", 
                choices = c('Leukemia (ALL): BCR/ABL mutated vs wild type'='OurData', namesExampleFile()))
    # choices = c('Leukemia (ALL): BCR/ABL mutated vs wild type'='OurData', 
    #             "Alzheimer's Disease" = "GSE1297", 
    #             "Renal Cancer" = "GSE14762", 
    #             "Pancreatic Cancer" = "GSE15471", 
    #             "Pancreatic Cancer" = "GSE16515", 
    #             "Non Small Cell Lung Cancer" = "GSE18842",
    #             "Glioma" = "GSE19728",
    #             "Colorectal cancer" = "GSE23878",
    #             "Thyroid Cancer" = "GSE3467", 
    #             "Alzheimer's Disease" = "GSE5281_EC",
    #             "Endometrial cancer" = "GSE7305", 
    #             "Acute myeloid leukemia" = "GSE9476"))
  })
  
  ## Download our example data set, loaded from sanssouci.data
  
  ### load data sets in the button : save in three csv file containing in a zip file
  output$downloadExampleData <- downloadHandler(
    filename = function() {
      paste("ExampleData", "zip", sep=".")
    },
    content = function(fname) {
      exampleData <- exampleData()
      fs <- c()
      tmpdir <- tempdir()
      setwd(tempdir())
      paths <- c("expressData.csv", "biologicalFunction.csv", "volcanoData.csv")
      write.csv(exampleData$matrix, paths[1])
      write.csv(exampleData$biologicalFunc, paths[2])
      write.csv(exampleData$degrade, paths[3])
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
    data(NULL)
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
    data(NULL)
    reset("fileGroup")
  })
  
  ## Cleaning data
  object_I <- 
    reactive({
      # eventReactive(input$buttonValidate, {
      withProgress(value = 0, message = "Upload Data... ", {
        if (input$checkboxDemo){ #if example data set
          if(req(input$choiceGSEA)=='OurData'){ # cleaning for data from sanssouci.data
            setProgress(value = 0.4, detail = "sanssouci data set ...")
            matrix <- expr_ALL #read data from sanssouci.data
            
            ### cleaning categories 
            cat <- colnames(matrix)
            categ <- rep(1, length(cat))
            categ[which(cat == "NEG")] <- 0
            
            object <- SansSouci(Y = as.matrix(matrix), groups = as.numeric(categ)) #create SansSouci object
            object$input$geneNames <- rownames(matrix)
            
            setProgress(value = 0.7, detail = "Preparation of gene set data ...  ")
            
            bioFun <- expr_ALL_GO
            stopifnot(nrow(bioFun) == nrow(matrix))  ## sanity check: dimensions
            mm <- match(base::rownames(bioFun), base::rownames(matrix))
            stopifnot(!any(is.na(mm)))
            object$input$biologicalFunc <- bioFun[mm, ]
            # print(dim(object$input$biologicalFunc)[1])
            rm(bioFun)
            object$bool$validation <- TRUE 
            object$bool$degrade <- FALSE
            rm(matrix)
            rm(categ)
          } else { # cleaning data set from GSEA data set
            setProgress(value = 0.4, detail = "GSEA data set ...")
            rawData <- readRDS(paste("GSEABenchmarkeR/express-data-set/", input$choiceGSEA, ".RDS", sep=""))
            # rawData <- R.cache::memoizedCall(maPreproc,geo2kegg()[input$choiceGSEA])[[1]]
            
            matrix <- SummarizedExperiment::assays(rawData)$exprs
            
            cats <- SummarizedExperiment::colData(rawData)
            ww <- match(cats$Sample, base::colnames(matrix))
            categ <- cats$GROUP[ww]
            object <- SansSouci(Y = matrix, groups = as.numeric(categ))
            setProgress(value = 0.7, detail = "GSEA data set ...")
            
            object$input$geneNames <- base::rownames(matrix)
            
            object$input$biologicalFunc <- readRDS("GSEABenchmarkeR/gene-set/go.gs.RDS")
            #On a laissé sous forme de liste car on a adapté les fonctions qui en ont besoin. Plus rapide qu'en la transformant en matrice binaire
            # print(length(object$input$biologicalFunc))
            object$bool$url <- rawData@metadata$experimentData@url
            
            object$bool$validation <- TRUE
            object$bool$degrade <- FALSE
          }
          
        } else { # if user use his own data
          req(input$fileData)
          req(fileData())
          file <- req(fileData())
          setProgress(value = 0.1, detail = "Read csv ...")
          matrix <- readCSV_sep(file = file$datapath, row.names = 1, check.names=FALSE)
          
          boolDegrade <- (length(matrix) == 2 & all(sort(colnames(matrix)) == sort(c("fc","p.value")))) #est ce que la matrice est dégradée
          setProgress(value = 0.2, detail = "Test data ...")
          if (boolDegrade){ #cleaning degraded matrix
            setProgress(value = 0.4)
            
            m=dim(matrix)[1]
            input <- list(m = m, geneNames = base::rownames(matrix))
            alpha = 0.1
            parameters <- list(
              B = 0,
              family = 'Simes', k=m)
            output <- list(p.value = matrix[["p.value"]], estimate = matrix[["fc"]])
            
            object <- structure(list(input = input,
                                     parameters = parameters,
                                     output = output), 
                                class = "SansSouci")
            
            object$bool$validation <- TRUE
            
          } else { #cleaning expression matrix
            setProgress(value = 0.4, detail = "Clean data set ...")
            clean <- cleanMatrix(matrix) #see cleanMatrix function : 
            ### gve specific issue for non available matrix
            
            ### if matrix is not available, matrix == NULL allows to block the calibration JER and the printing  
            if(clean$boolValidation){ #if matrix is ok
              object <- SansSouci(Y = as.matrix(clean$data), groups = as.numeric(colnames(clean$data)))
              object$input$geneNames <- base::rownames(clean$data)
            }else{ #if isn't 
              input <- NULL
              parameters <- NULL
              output <- NULL
              
              object <- structure(list(input = input,
                                       parameters = parameters,
                                       output = output), 
                                  class = "SansSouci")
            }
            
            object$bool$validation <- clean$boolValidation
            object$bool$matrix.color <- clean$color #color of error message
            object$bool$matrix.text <- clean$text #eroor message
            
          }
          
          object$bool$degrade <- boolDegrade
          
          setProgress(value = 0.7, detail = "Preparation of gene set data ...")
          
          ## cleaning of gene set matrix
          fileGroup <- fileGroup()
          if (!is.null(fileGroup)){
            setProgress(value = 0.75, detail = "Read gene set data ...")
            T1 <- Sys.time()
            bioFun <- readCSV_sep(file = fileGroup$datapath, row.names = 1, check.names=FALSE)
            T2 <- Sys.time()
            print(paste("read gene set", T2-T1))
            setProgress(value = 0.8, detail = "Cleaning of gene set data ...")
            cleanBio <- cleanBiofun(bioFun) # cleaning and message error
            
            rm(bioFun)
            
            object$bool$bioFun.color <- cleanBio$color
            object$bool$bioFun.text <- cleanBio$text
            
            
            setProgress(value = 0.9, detail = "Matching ...")
            
            matchBio <- matchMatrixBiofun(geneNames = object$input$geneNames, biofun = cleanBio$biofun) #verification compatibility between the two matrices
            if(matchBio$boolValidation & cleanBio$boolValidation){ #if ok
              object$input$biologicalFunc <- as.matrix(matchBio$biofun)
            } else { #if not ok
              object$input$biologicalFunc <- NULL
            }
            object$bool$match.color <- matchBio$color
            object$bool$match.text <- matchBio$text
            
            rm(matchBio)
            rm(cleanBio)
            
          }
          
          rm(matrix)
          
          
          
        }
        setProgress(value = 1, detail = "Done")
        return(object)
      })
    }
    )
  
  data <- reactiveVal()
  observe({
    data(req(object_I()))
  })
  
  
  
  urlDataSet <- eventReactive(input$buttonValidate, { #give link to description of geo2kegg data sets
    # req(geo2kegg())
    req(input$choiceGSEA)
    # req(geo2kegg()[[input$choiceGSEA]])
    req(object_I()$bool$url)
    # return(a("URL link to data set description", href=geo2kegg()[[input$choiceGSEA]]@experimentData@url))
    return(a("URL link to data set description", href=object_I()$bool$url, target="_blank"))
    
  })
  output$msgURLds <- renderUI({ 
    req(urlDataSet)
    tagList(urlDataSet())
  })
  
  ## different error message for non compliant matrix
  output$errorInput <- renderUI({ # express gene matrix error
    tags$span(style= req(object_I()$bool$matrix.color), paste(req(object_I()$bool$matrix.text)))
  })
  
  output$errorBioMatrix <- renderUI({ #gene set matrix error
    tags$span(style=req(object_I()$bool$bioFun.color), paste(req(object_I()$bool$bioFun.text)))
  })
  output$errorMatch <- renderUI({ # none commun gene btw expression gene matrix and gene set matrix
    tags$span(style=req(object_I()$bool$match.color), paste(req(object_I()$bool$match.text)))
  })
  
  
  output$watch <- renderPrint({tableCSV()})
  
  
  
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
    req(object_I())
    numericInput("valueK", 
                 label = "K (size of reference family)", 
                 value = numKI(),
                 min = 1,
                 max = nrow(req(object_I()$input$Y)))
  })
  
  ## numKI() is used to intiate the printed input valueK 
  numKI <- reactiveVal()
  observe({ #Initialization, if refFamily == 'Beta' or not 
    req(object_I())
    newValue <- ifelse(input$refFamily == "Beta", 
                       round(2*req(nrow(object_I()$input$Y))/100),
                       req(nrow(object_I()$input$Y)))
    numKI(newValue)
  })
  
  ## numK is the parameters choosen by users and use in server side
  numK <- reactiveVal()
  observeEvent(object_I(),{   #si les paramètres ne sont pas activés, input$valueK n'existe pas et donc on a pas le bon résultat ... pourquoi ? CalibrateJER devrait gérer si K est null non ? 
    req(object_I())           # when object_I() is change, INITIALISATION of numK()
    newValue <- req(nrow(object_I()$input$Y))
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
                                         "Reference family = 'Simes',\n K = ", length(object_I()$output$p.value)))
  })
  
  observe({
    req(object_I())
    if(object_I()$bool$degrade){
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
  
  
  
  # ## JER calibration on available expression gene matrix
  # cal <- reactive({
  #   # eventReactive(input$buttonValidate, {
  #   withProgress(value = 0, message = "Perform calibration ... ", {
  #     incProgress(amount = 0.3)
  #     t1 <- Sys.time()
  # 
  #     cal <- R.cache::memoizedCall(calibrateJER,
  #                                  req(data()$matrix), # if data()$matrix == NULL, not perform [-> not available matrix or degraded matrix]
  #                                  categ = data()$categ,
  #                                  B = numB(), alpha = alpha(),
  #                                  refFamily = refFamily(), alternative = alternative(),
  #                                  K = numK()
  #     )
  #     
  #     
  #     t2 <- Sys.time()
  #     print(paste("calibration :",difftime(t2, t1)))
  #     setProgress(value = 0.7, detail = "Done")
  #     return(cal)
  #   })
  # })
  ## JER calibration on available expression gene matrix
  # observe({
  observeEvent(input$buttonValidate, {
    req(data()$bool$validation)
    if(!data()$bool$degrade){ #non degraded version
      withProgress(value = 0, message = "Perform calibration ... ", {
        incProgress(amount = 0.3)
        t1 <- Sys.time()
        # object <- memoizedCall(fit, 
        #                        object = data(), 
        #                        alpha = req(alpha()),
        #                        B = numB(),
        #                        alternative = alternative(),
        #                        family = refFamily(),
        #                        K = numK()
        # )
        object <- fit(data(),
                      alpha = req(alpha()),
                      B = numB(),
                      alternative = alternative(),
                      family = refFamily(),
                      K = numK()
        )
        
        
        t2 <- Sys.time()
        print(paste("calibration :",difftime(t2, t1)))
        setProgress(value = 0.7, detail = "Done")
        # data(object)
      })
    } else { #degraded version #if degraded matrix is used
      #les matrices pval et fc sont déjà rentrées dans l'objet lors que la création de l'objet
      
      # object <- fit(data(),   #Non utilisable car nécessite Y
      #               alpha = req(alpha()),
      #               B = 0,
      #               alternative = alternative(),
      #               family = refFamily(),
      #               K = numK()
      # )
      
      m = nHyp(data())
      thr <- t_linear(alpha(), seq_len(m), m) # force using of Simes and k=m # IMPORT FROM FUNCTION.R
      object <- data()
      object$output$thr <- thr
      object$output$lambda <- alpha()
      object$parameters$alpha <- alpha()
      # data(object)
    }
    #calcul des logp et adjp
    object$output$logp <- -log10(pValues(object))
    object$output$adjp <- p.adjust(pValues(object), method = "BH")
    data(object)
  })
  
  ###### dire ce qui est calculé dans cette partie ci dessus
  
  # ici sont calculés la calibration (thr et stat pivotale) 
  # puis rowWelchTest : -> obtient p-values ('p.value') et foldchange ('estimate')
  
  
  
  ## Thr calculation
  # thr <- reactiveVal() #Initialization
  # 
  # observe({ 
  #   req(data()$matrix) # if gene expression matrix is available 
  #   newValue <- req(cal()$thr) #we use thr from calibration and user's parameters
  #   thr(newValue)
  # })
  # 
  # observe({ # if p.value matrix is used
  #   
  #   
  #   
  #   req(data()$df) #if degraded matrix is used
  #   # if the gene expression matrix is available => data()$df == NULL 
  #   m = dim(data()$df)[1]
  #   newValue <- SimesThresholdFamily(m, kMax = m)(alpha()) # force using of Simes and k=m
  #   thr(newValue)
  # })
  
  ################### 
  # P-values and foldchange
  ###################
  
  # df <- reactiveVal() #creation
  
  # observe({ # If gene expression matrix is available
  #   req(data())
  #   req(pValues(data()))
  #   req(thresholds(data()))
  #   object <- req(data()) # on récupère l'objet
  #   object$output$logp <- -log10(pValues(object))
  #   object$output$adjp <- p.adjust(pValues(object), method = "BH")
  #   data(object)
  #   # req(data()$matrix) #active observe cell
  #   # dex <- rowWelchTests(req(data()$matrix), data()$categ) # make a test
  #   # print(dex)
  #   # pval <- dex[["p.value"]]
  #   # logp <- -log10(pval)
  #   # fc <- dex$estimate
  #   # adjp <- p.adjust(pval, method = "BH")
  #   # newValue <- list(logp = logp, fc = fc, adjp = adjp, pval = pval)
  #   # df(newValue) #update df()
  # })
  
  # observe({ # if degraded version
  #   req(data()$df) #active observe cell
  #   pval <- data()$df[['p.value']] # from input
  #   fc <- data()$df[['fc']]
  #   logp <- -log10(pval)
  #   adjp <- p.adjust(pval, method = "BH")
  #   newValue <- list(logp = logp, fc = fc, adjp = adjp, pval = pval)
  #   df(newValue)
  # })
  
  
  ################### 
  # Threshold 
  ###################
  
  # vertical contient la valeur de déplacement de tous les objets déplaçable (dans notre cas juste les seuils)
  vertical <- reactive({event_data("plotly_relayout", source='A')})
  
  ###
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
  observeEvent(data()$output$p.value, { # initialisation : adjusted p-value == 0.1
    req(data())
    min0.1 <- which.min(abs(data()$output$adjp-0.1))
    y0.1 <- data()$output$logp[[min0.1]]
    # print(paste("min0.1 = ",min0.1," y0.1 = ", y0.1))
    yint(y0.1)
  })
  observeEvent(vertical()[["shapes[1].y0"]], { # when user change threshold on plotly
    p <- 10^(-vertical()[["shapes[1].y0"]])
    y_sel <- which((pValues(data()) <= p))          ## selected by  p-value
    newValue <- Inf
    if (length(y_sel) > 0) {
      newValue <- min(data()$output$logp[y_sel])              ## threshold on the log(p-value) scale
    }
    yint(newValue)
  })
  
  ################### 
  # Selecting gene with thresholds
  ###################
  
  # selected genes for server calcuation
  selectedGenes <- reactive({
    req(data()$output$logp)
    ## gene selections
    sel1 <- which(data()$output$logp >= yint() & data()$output$estimate >= xint()) #upper right 
    sel2 <- which(data()$output$logp >= yint() & data()$output$estimate <= xint2()) #upper left 
    sel12 <- sort(union(sel1,sel2)) #both
    return(list(sel1 = sel1, sel2 = sel2, sel12 = sel12))
  })
  
  
  #matrix p-val & fc of selected genes by thresholds => to plot red points on volcano plot
  selected_points <- reactive({
    list(x = data()$output$estimate[selectedGenes()$sel12], 
         y = data()$output$logp[selectedGenes()$sel12])
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
    ###
    ### both
    # n12 <- length(selectedGenes()$sel12)
    # FP12 <- maxFP(df()$pval[selectedGenes()$sel12], thr = thr())
    # TP12 <- n12 - FP12
    # FDP12 <- round(FP12/max(n12, 1), 2)
    
    req(thresholds(data()))
    n12 <- length(selectedGenes()$sel12)
    pred <- predict(object = data(), S = selectedGenes()$sel12, what = c("TP", "FDP"))
    return(list(
      n12 = n12, TP12 = pred["TP"], FDP12 = pred["FDP"]))
    
  })
  
  # calculate PHB for lasso box selection
  calcBoundSelection <- reactive({ # 
    req(manuelSelected())
    req(thresholds(data()))
    c(n = length(manuelSelected()), predict(data(), S=manuelSelected(), what = c("TP", "FDP")))
    # calcBounds(df()$pval[manuelSelected()], thr = thr()) #return(list(n=n, FP=FP, TP=TP, FDP=FDP))
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
               "FDP≤" = c(round(TP_FDP()$FDP12, 2)),
               check.names = FALSE, row.names = NULL)
  })
  
  # updating PHB table when thresholds change
  observeEvent(TP_FDP(),{ # When threshold change
    
    
    bottomTable <- tableResult() %>%  #keep box/lasso selection
      filter(Selection != "Threshold selection") #remove row named "Threshold selection"
    upperTable <- baseTable() # take new value of PHB for the new threshold selection
    newValue <- rbind(upperTable, bottomTable)
    tableResult(newValue)
  })
  
  # Add PHB for lasso and box selection 
  ## d() is the reactive variable for box and lasso selection
  observeEvent(d(),{  # When user selects a new group of points
    req(calcBoundSelection())
    vectorGene <- names(pValues(data())[manuelSelected()]) #list of gene contained in gene selection
    url <- UrlStringdbGrah(vectorGene) #construction of link to StringDB interaction graph
    n <- dim(tableResult())[1]
    newValue <- rbind(tableResult(), c(paste('<a target="_blank" href="', url, '" >User selection ',n, '</a>', sep=""), 
                                       calcBoundSelection()['n'], 
                                       calcBoundSelection()['TP'],
                                       round(calcBoundSelection()['FDP'], 2)))
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
    req(data()$output$adjp)
    listLog <- c()
    for (i in c(0.5,0.25,0.1,0.05,0.025,0.01,0.001,0.0001)){ #selected line on yaxis /!\ if you change it you have to change in yaxis()
      min05 <- which.min(abs(data()$output$adjp-i)) #to take the gene with the nearest adjpvalue to i 
      y05 <- data()$output$logp[[min05]] # take the equivalent in logpvalue to print it on VP
      listLog <- c(listLog, y05)
    }
    return(listLog)
  })
  
  # value for 'NUmber of false discoverie' yaxis : optimize y line
  thr_yaxis <- reactive({
    req(alpha()) # if parameters change (not only alpha)
    req(data())
    req(thresholds(data()))
    req(data()$output$logp)
    thrYaxis(thr = thresholds(data()), maxlogp=max(data()$output$logp))
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
    req(data()$bool$validation)
    
    req(foldChanges(data()))
    # req(data()$output$logp)
    # print('posteriori()')
    setProgress(value = 0.9, detail = "posteriori Reactive ... ")
    f <- list(
      size = 14,
      color = "#000000"
    )
    lte <- "≤"
    gte <- "≥"
    plot_ly(data.frame(x = foldChanges(data()), y=data()$output$logp), 
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
    vect <- tableCSV()[,req(userDTselectPost())] 
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
            x = unname(foldChanges(data())[selectionUserRe()$sel]),
            y = unname(data()$output$logp[selectionUserRe()$sel]),
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
  observeEvent({
    input$resetCSV
    data()$input$Y
    input$buttonValidate
  },{
    req(data()$input$Y)
    req(selectedGenes())
    vecteur <- rep(0, dim(data()$input$Y)[1])
    vecteur[selectedGenes()$sel12] <- 1
    newValue <- data.frame(Thresholds_selection = vecteur, row.names = req(data()$input$geneNames)) #a dataframe with rowname without features
    tableCSV(newValue)
  })
  
  #When user want to change thresholds
  observeEvent(selectedGenes(),{
    req(data()$input$Y)
    req(selectedGenes())
    req(tableCSV())
    vecteur <- rep(0, dim(data()$input$Y)[1])
    vecteur[selectedGenes()$sel12] <- 1
    
    rigthTable <- tableCSV() %>% 
      select(-Thresholds_selection)
    df <- cbind(data.frame(Thresholds_selection = vecteur, row.names = req(data()$input$geneNames)), 
                rigthTable)
    tableCSV(df)
  })
  
  #When user selected ne gene set from lasso/box [d() activate]
  observeEvent(d(), { 
    req(data()$input$Y)
    req(tableCSV())
    vecteur <- rep(0, dim(data()$input$Y)[1])
    vecteur[manuelSelected()] <- 1
    nameCol <- colnames(tableCSV())
    df <- cbind(tableCSV(), selection2=vecteur)
    colnames(df) <- c(nameCol, paste("User selection", length(nameCol)))
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
      table <- tableResult()
      table$`Selection` <- str_remove_all(str_remove_all(table$`Selection`, "<a(.*?)>"), "(</a>)")
      write.csv(table, file)
    }
  )
  
  # while buttonValidate is not clicke, these button are hidden
  observeEvent(input$buttonValidate, {
    shinyjs::show("downloadPHBTable")
    shinyjs::show("resetCSV")
    shinyjs::show("downloadData")
  })
  
  # /!\ si on réactive ces outputs il faut les mettre à jour avec la nouvelle version OBJECT de SANSSOUCI
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
  
  # output$outThresholds <- renderPrint({
  #   print("Thresholds(data()")
  #   
  #   print(req(thresholds(req(data())))[1])
  # })
  # output$outbioFun <- renderPrint({
  #   print("data()$input$biologicalFunc")
  #   
  #   print(req(data()$input$biologicalFunc)[1])
  # })
  # 
  # output$outThrBioFun <- renderPrint({
  #   print("data()$input$biologicalFunc")
  #   req(thresholds(req(data())), data()$input$biologicalFunc)
  #   print(req(thresholds(req(data())))[1])
  #   print(req(data()$input$biologicalFunc)[1])
  # })
  
  # PHB for all gene sets
  tableBoundsGroup <- reactive({
    withProgress(message = "tableBoundsGroup", {
      T1 <- Sys.time()
      req(thresholds(req(data())))
      req(data()$input$biologicalFunc)
      # print('tableBoundsGroup()')
      # print(str(data()))
      # req(data()$input$biologicalFunc)
      table <- boundGroup2(req(data())
      )
      T2 <- Sys.time()
      # print(paste("post hoc bounds on gene set:",difftime(T2, T1)))
      print("post hoc bounds on gene set:")
      print(difftime(T2, T1))
    })
    return(table)
  })
  
  # calculate vounds for all features to compare with each gene sets (competitive methods)
  boundsW <- reactive({ #calcul bounds for all features
    req(pValues(data()))
    req(thresholds(data()))
    c(n=length(nHyp(data())), predict(object = data()))
  })
  
  
  # dowload csv file containing PHB table of gene sets
  output$downloadPHBTableGroup <- downloadHandler( #download csv of user selection
    filename = function() {
      tag <- format(Sys.time(), "%Y-%M-%d_%H-%m-%S")
      sprintf("gene-set_bounds_%s.csv", tag)
    },
    content = function(file) {
      table <- tableBoundsGroup()
      table$Name <- str_remove_all(str_remove_all(table$Name, "<a(.*?)>"), "(</a>)")
      write.csv(table, file)
    }
  )
  
  # observe({ # show the download button when 
  #   req(tableBoundsGroup())
  #   shinyjs::show("downloadPHBTableGroup")
  # })
  
  # If user choose SEA alternative
  filteredTableBoundsGroup  <- reactive({
    req(input$buttonSEA)
    req(tableBoundsGroup())
    if (input$buttonSEA == "competitive"){
      table <- tableBoundsGroup()
      sel <- which(table[["FDP≤"]] < boundsW()['FDP'])
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
           title = "Data", content = paste("This table prints your post-hoc bounds for your gene sets."  ,
                                           "For example, the selection called",
                                           tableBoundsGroup()[1, "Name"],
                                           "contains at leat", 
                                           tableBoundsGroup()[1, "TP≥"],
                                           " true positives (TP) and its False Discovery Proportion (FDP) is less than ",
                                           round(tableBoundsGroup()[1, "FDP≤"]*100, 2), "%"
           ), 
           trigger='hover')
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
    # req(input$tableBoundsGroup_rows_selected)
    href <- filteredTableBoundsGroup()[input$tableBoundsGroup_rows_selected,"Name"]
    name <- str_remove_all(str_remove_all(href, "<a(.*?)>"), "(</a>)")
    return(name)
  })
  
  #list of genes selected by user
  selectionGroup <- reactive({
    req(data()$input$biologicalFunc)
    
    group <- req(userDTselectPrio())
    bioFun <- data()$input$biologicalFunc
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
    req(data())
    req(foldChanges(data()))
    req(data()$output$logp)
    f <- list(
      size = 14,
      color = "#000000"
    )
    lte <- "≤"
    gte <- "≥"
    print("Priori()")
    plot_ly(data.frame(x = foldChanges(data()), y=data()$output$logp), 
            x = ~x, y = ~y, 
            marker = list(size = 2,
                          color = 'grey'), 
            name = 'genes',
            type='scattergl', mode = "markers", source='B'
            ,
            text = data()$input$geneNames,
            customdata = paste0("http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=", data()$input$geneNames)
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
            x = unname(foldChanges(data())[selectionGroup()$sel]),
            y = unname(data()$output$logp[selectionGroup()$sel]),
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
