library("shiny")
library("plotly")
m=1000
xi=rnorm(m)
#theta<-rep(0, m)
#pbern<-3/4
# a=0.95
# b=0.65
# A<-matrix(c(a, 1-a, 1-b, b), 2, 2, byrow=T)
# theta[1]=0
# for (i in 2:m)
# { theta[i] <- (1-theta[i-1])*rbinom(1, 1, A[1, 2]) + theta[i-1]*rbinom(1, 1, A[2, 2])
# }
ChooseWi<-"100px"
ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            #numericInput("scale", "Resolution : ", 64, min = 32, max = 128),
            numericInput("alpha", "Target confidence level:", 0.1, min = 0, max = 1,width=ChooseWi),
            numericInput("delta", "Signal strengh:", 3, min =1 , max = 10,width=ChooseWi),
            numericInput("s", "Region size:", 10, min =1 , max = m/10,width=ChooseWi),
            numericInput("nbRegsignal", "Number of regions with signal:", 3, min =1 , max = 10,width=ChooseWi),
            #numericInput("sss", "Random seed:",1,min=1,max=100000,width=ChooseWi),
            width = 2
        ),
        mainPanel(
            h3(textOutput("bound")),
            plotlyOutput("plot")
        )
    )
)
server <- function(input, output, session) {
    
    re <- reactive({
        moy=input$delta
        s=input$s
        nbRegsignal=input$nbRegsignal
        #set.seed(input$sss)
        #s=200
        K=m/s
        #nbRegsignal=3
        Regsignal=sample(K,nbRegsignal)
        signal=as.numeric(sapply(Regsignal,function(x){(1:s)+(x-1)*s}))
        theta=rep(0,m)
        theta[signal]=1
        #theta<-theta*rbinom(m,1,pbern)
        muvector=theta*moy
        X=muvector+xi
        pvalues=pnorm(X,lower.tail=FALSE)
        list(pvalues = pvalues, Regsignal = Regsignal, K = K)
    })
    
    output$plot <- renderPlotly({
        res <- re()
        pvalues <- res$pvalues
        
        data=data.frame(key=1:m,SNP=1:m,mlogp=-log(pvalues))
        p=plot_ly(data,  x=~SNP, y = ~mlogp, key = ~key,width=800,height=400) %>%
            layout(dragmode = "select")  %>% layout(xaxis= list(title=""),yaxis= list(title=""))
        #      add_markers(p, color = ~signal)
        
    })
    
    output$bound <- renderText({
        s=input$s
        alpha=input$alpha
        res <- re()
        pvalues <- res$pvalues
        Regsignal <- res$Regsignal
        K <- res$K
        
        d <- event_data("plotly_selected")
        if (is.null(d)) {
            collec=c(Regsignal)
            R=unique(as.numeric(sapply(collec,function(x){(1:s)+(x-1)*s})))
        }else{
            R=as.numeric(d$key)
        }
        pval=pvalues[R]
        
        V=min(sapply(1:m,function(k) min(sum(pval>=alpha*k/m)+k-1,length(pval) ) )) #/length(pval)
        # truc=function(k){
        #     Rk=((k-1)*s+1):(k*s)
        #     #zetak=sum(pvalues[Rk]>alpha/m)
        #     zetak=s-HB(pvalues[Rk],alpha/K)$hatk
        #     return(min(zetak,length(intersect(Rk,R))))
        # }
        #Vpart=sum(sapply(1:K,truc))  #/length(pval)
        truc2=function(k){
            Rk=((k-1)*s+1):(k*s)
            ualpha=sqrt(log(K/alpha)/2)
            zetak=sapply(pvalues,function(t) 
                (sqrt(sum(pvalues[Rk]>t)/(1-t)+(ualpha/(2*(1-t)))^2) + ualpha/(2*(1-t)))^2)
            #zetak=min(sapply(pvalues,function(t)
            # sum(pvalues[Rk]>t)/max(10^(-16),1-t-sqrt(log(K/alpha)/(2*length(Rk))))))
            return(min(zetak,length(intersect(Rk,R))))
        }
        VDKW = floor(sum(sapply(1:K,truc2)))
        
        
        # if (is.null(d)) {
        #   msg <- "Select a set of points"
        # } else {
          #msg <- sprintf("FDP smaller than %s or %s among %s selected pixels", signif(V/length(pval),4), signif(VDKW/length(pval),4), length(R))
          msg <- sprintf("True nulls smaller than %s (Simes) or %s (part) among %s selected pixels", V, VDKW, length(R))  
        # }
        msg
    })
}
shinyApp(ui, server)




