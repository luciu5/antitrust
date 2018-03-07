require(shiny)
require(rhandsontable)


shinyServer(function(input, output, session) {


   nPossProds <- 10 #only allow 10 products
   
   msgCatcher <- 

   function(expr)
   {
     W <- NULL
     E <- NULL
     
     w.handler <- function(w){ # warning handler
       W <<- append(W,conditionMessage(w))
       #invokeRestart("muffleWarning")
     }
     e.handler <- function(e){ # error handler
       E <<- append(E, conditionMessage(e))
       NULL
       
     }
     list(value = withCallingHandlers(tryCatch(expr, error = e.handler),
                                      warning = w.handler),
          warning = W, error = E)
   }
   
   
   
   
   genInputData <- function(){
    # a function to generate default input data set for simulations
     
     inputData <- data.frame(
       Name = c("Prod1","Prod2","Prod3","Prod4"),
       ownerPre  = c("Firm1","Firm2","Firm3","Firm3"),
       ownerPost = c("Firm1","Firm1","Firm3","Firm3"),
       'Prices \n($/unit)'    =c(.0441,.0328,.0409,.0396)*100,
       Shares   =c(0.09734513, 0.25368732, 0.37315634, 0.27581121),
       Margins =c(.3830,.5515,.5421,.5557),
       stringsAsFactors = FALSE,
       check.names=FALSE
     )
     
     
     
     nDefProd <- nrow(inputData)
     inputData <- inputData[c(1:nDefProd,rep(1, nPossProds - nDefProd)),]
     #if(input$incEff) inputData$mcDelta <- 0
     
     inputData[(nDefProd+1):nPossProds,] <- NA
     rownames(inputData) <- NULL
     
    
     return(inputData)

   }
   
   gensum <- function(res){
     #a function to generate summary stats for results tab
     
     
       
     isCournot <- grepl("Cournot",class(res))
     isRevDemand <- grepl("ces|aids",class(res),ignore.case = TRUE)
     
     if(isCournot){ 
      
       capture.output( s <- summary(res, market = FALSE))
        s$sharesPost <- s$quantityPost/sum(s$quantityPost,na.rm=TRUE)*100
     }
     
     else{capture.output(s <- summary(res))}
     
     thiscmcr <- thisCV <- NA
     try(thiscmcr <- cmcr(res), silent=TRUE)
     try(thiscv <- CV(res),silent = TRUE)
     #thiselast <- elast(res,market=TRUE)
     
     
     thispsdelta  <- sum(calcProducerSurplus(res,preMerger=FALSE) -calcProducerSurplus(res,preMerger=TRUE),na.rm=TRUE)
     
     ## $ is the symbol for generic currency. replace in favor of $
     res <- data.frame(
           'HHI Change' = as.integer(round(hhi(res,preMerger=FALSE) -  hhi(res,preMerger=TRUE))),
           'Industry Price Change (%)' = sum((s$priceDelta)*s$sharesPost/100,na.rm=TRUE),
           'Merging Party Price Change (%)'= sum(((s$priceDelta)*s$sharesPost)[s$isParty == "*"] / s$sharesPost[s$isParty== "*"],na.rm=TRUE),
           'Compensating Marginal Cost Reduction (%)' = ifelse(isCournot, thiscmcr, sum(thiscmcr*s$sharesPost[s$isParty== "*"]/sum(s$sharesPost[s$isParty== "*"],na.rm=TRUE))),
           'Consumer Harm ($/unit)' = -1*thiscv,
           'Producer Benefit ($/unit)' = thispsdelta,
           'Overall Effect ($/unit)'= -1*thiscv + thispsdelta,
           
           #'Estimated Market Elasticity' = thiselast,
           check.names=FALSE
     )
     
     if(isRevDemand){
       res$'Overall Effect ($/unit)' <- NULL
       colnames(res) <- gsub('(?<=Consumer Harm\\s)\\(\\$/unit\\)',"(% Expenditure)",colnames(res), perl=TRUE)
       }
     
     
     if(is.na(res[,"Compensating Marginal Cost Reduction (%)"])) res[,"Compensating Marginal Cost Reduction (%)"] <- NULL
     
     return(res)
   }


   
   runSims <- function(supply, demand, indata, mktElast){
     # a function to execute code from antitrust package based on ui inputs
     
     firstMargin <- which(!is.na(indata$Margins))[1]
     firstPrice <- which(!is.na(indata[,"Prices \n($/unit)"]))[1]
     
     
     ownerPre = model.matrix(~-1+indata$ownerPre)
     ownerPre = tcrossprod(ownerPre)
     ownerPost = model.matrix(~-1+indata$ownerPost)
     ownerPost = tcrossprod(ownerPost)
     
     
     switch(supply,
            Bertrand =
              switch(demand,
                     `logit (unknown elasticity)`= logit.alm(prices= indata[,"Prices \n($/unit)"],
                                                             shares= indata[,"Output"],
                                                             margins= indata$Margins,
                                                             ownerPre= ownerPre,
                                                             ownerPost= ownerPost,
                                                             mcDelta = indata$mcDelta, labels=indata$Name),
                     `ces (unknown elasticity)`= ces.alm(prices= indata[,"Prices \n($/unit)"],
                                                         shares= indata[,"Output"],
                                                         margins= indata$Margins,
                                                         ownerPre= ownerPre,
                                                         ownerPost= ownerPost,
                                                         mcDelta = indata$mcDelta, labels=indata$Name),
                     linear=linear(prices= indata[,"Prices \n($/unit)"],
                                   quantities= indata[,"Output"],
                                   margins= indata$Margins,
                                   ownerPre= ownerPre,
                                   ownerPost= ownerPost,
                                   mcDelta = indata$mcDelta, labels=indata$Name),
                     aids=aids(prices= indata[,"Prices \n($/unit)"],
                               shares= indata[,"Output"],
                               margins= indata$Margins,
                               ownerPre= ownerPre,
                               ownerPost= ownerPost,
                               mcDelta = indata$mcDelta, labels=indata$Name),
                     logit= logit.alm(prices= indata[,"Prices \n($/unit)"],
                                      shares= indata[,"Output"],
                                      margins= indata$Margins,
                                      ownerPre= ownerPre,
                                      ownerPost= ownerPost,
                                      mcDelta = indata$mcDelta, labels=indata$Name,  mktElast = mktElast ),
                     ces = ces.alm(prices= indata[,"Prices \n($/unit)"],
                                   shares= indata[,"Output"],
                                   margins= indata$Margins,
                                   ownerPre= ownerPre,
                                   ownerPost= ownerPost,
                                   mcDelta = indata$mcDelta, labels=indata$Name,  mktElast = mktElast),
                     linear=linear(prices= indata[,"Prices \n($/unit)"],
                                   quantities= indata[,"Output"],
                                   margins= indata$Margins,
                                   ownerPre= ownerPre,
                                   ownerPost= ownerPost,
                                   mcDelta = indata$mcDelta, labels=indata$Name),
                     pcaids=pcaids(prices= indata[,"Prices \n($/unit)"],
                                   shares= indata[,"Output"],
                                   knownElast = -1/indata$Margins[firstMargin],
                                   knownElastIndex = firstMargin,
                                   mktElast = mktElast,
                                   ownerPre= ownerPre,
                                   ownerPost= ownerPost,
                                   mcDelta = indata$mcDelta, labels=indata$Name)
              ),
            Cournot = 
              
              cournot(prices= indata[,"Prices \n($/unit)"][firstPrice],
                      demand = gsub("\\s+\\(.*","",demand,perl=TRUE),
                      cost= rep("linear", nrow(indata)),
                      quantities = as.matrix(indata[,"Output"]),
                      margins= as.matrix(indata$Margins),
                      ownerPre= ownerPre,
                      ownerPost= ownerPost,
                      mktElast = ifelse( grepl("unknown elasticity", demand),
                                         NA_real_, mktElast),
                      mcDelta = indata$mcDelta, 
                      labels=list(as.character(indata$ownerPre),indata$Name[firstPrice]))
            ,
            `2nd Score Auction`= switch(demand,
                                        `logit (unknown elasticity)` = auction2nd.logit.alm(prices= indata[,"Prices \n($/unit)"],
                                                                                            shares= indata[,"Output"],
                                                                                            margins= indata$Margins,
                                                                                            ownerPre= ownerPre,
                                                                                            ownerPost= ownerPost,
                                                                                            mcDelta = indata$mcDelta, labels=indata$Name),
                                        logit = auction2nd.logit.alm(prices= indata[,"Prices \n($/unit)"],
                                                                     shares= indata[,"Output"],
                                                                     margins= indata$Margins,
                                                                     ownerPre= ownerPre,
                                                                     ownerPost= ownerPost,
                                                                     mcDelta = indata$mcDelta, labels=indata$Name, 
                                                                     mktElast = mktElast)
                                        
            )
     )
   }
   
   
   
   ## create a reactive object that tracks which demand is being used
   demand <- eventReactive(input$simulate, {ifelse(input$supply =="Bertrand", 

                                              ifelse(input$calcElast == "at least 2 margins",input$demand_bert_alm,input$demand_bert),
                                                      ifelse(input$supply == "Cournot",
                                                     ifelse(input$calcElast == "at least 2 margins",input$demand_cournot_alm, input$demand_cournot),
                                                     ifelse(input$calcElast == "at least 2 margins",input$demand_2nd_alm, input$demand_2nd)))
                      })
   
   
   ## create a reactive list of objects
   values <- reactiveValues(inputData = genInputData(), sim = NULL, msg = NULL)
      
   
                        
   ## update reactive list whenever changes are made to input
    observe({
    
      
      if(!is.null(input$hot)){
        values$inputData = hot_to_r(input$hot)
        
      }
    })
    

    ## display inputs 
    output$hot <- renderRHandsontable({
      
      inputData <- values[["inputData"]]
      
      if(input$supply =="Cournot"){colnames(inputData)[colnames(inputData) == "Shares"] <- "Quantities"}
      else{colnames(inputData)[colnames(inputData) == "Quantities"] <- "Shares"}
      
      if (!is.null(inputData))
        rhandsontable(inputData, stretchH = "all")
    })


    
  ## simulate merger when the "simulate" button is clicked
   observeEvent(input$simulate,{
      

     values[["sim"]] <- values[["msg"]] <-  NULL
     
     updateTabsetPanel(session,inputId  = "inTabset", selected = "respanel")
      
      indata <- values[["inputData"]]
      
      colnames(indata)[colnames(indata) %in% c("Quantities","Shares")] <- "Output"
      
      indata <- indata[!is.na(indata[,"Output"]),]
      
      
      
      indata$mcDelta <- 0
      
      indata$ownerPre <- factor(indata$ownerPre)
      indata$ownerPost <- factor(indata$ownerPost,levels=levels(indata$ownerPre))
      
      
     
      thisSim <- msgCatcher(
                          runSims(supply = input$supply,demand = demand(), indata = indata, mktElast = input$enterElast )
      )

     values[["sim"]] <-  thisSim$value
               
     values[["msg"]] <-  list(error=thisSim$error,warning=thisSim$warning)        
      
     if(!is.null(thisSim$error)) updateTabsetPanel(session,inputId  = "inTabset", selected = "msgpanel")
      
    })
    
    
   
      
    ## display summary results from gensum to results tab  
    output$results <-
          
          renderTable({
           
            if(input$inTabset != "respanel" || input$simulate == 0|| is.null(values[["sim"]])){return()}
         
           
            gensum(values[["sim"]])
          })
  
    
  ## display summary values to details tab      
    output$results_detailed <- renderTable({
       
      if(input$inTabset!= "detpanel" || input$simulate == 0  || is.null(values[["sim"]])){return()}
        
        if(input$supply == "Cournot"){
          
          res <- NULL
          capture.output(try(res <- summary(values[["sim"]], market=FALSE),silent=TRUE))
          
          
          res$product <- res$mcDelta <- NULL
          
          try(colnames(res) <- c("Party","Name","Pre-Merger Price","Post-Merger Price", "Price Change (%)","Pre-Merger Quantity","Post-Merger Quantity", "Output Change (%)"),silent=TRUE)
          
          }
        
        else{
          
          capture.output(res <- summary(values[["sim"]]))
          res$Name <- rownames(res)
          res <- res[,c(1, ncol(res), 2 : (ncol(res)  - 1))]
          colnames(res) <- c("Party","Name","Pre-Merger Price","Post-Merger Price", "Price Change (%)","Pre-Merger Share (%)","Post-Merger Share (%)", "Share Change (%)")
        }
        
        
        print(res)
        
      
      })  
    
    
   
    ## display elasticities to elasticity tab
    output$results_elast <- renderTable({
      
      if(input$inTabset!= "elastpanel" || input$simulate == 0 || is.null(values[["sim"]])){return()}
        
       if(input$pre_elast == "Pre-Merger"){ preMerger = TRUE}
        else{preMerger =FALSE}
        
        res <- elast(values[["sim"]], preMerger=preMerger)
        print(res)
        
      
    }, rownames = TRUE)   
  
    
    output$results_mktelast <- renderTable({
      
      if(input$inTabset!= "elastpanel" || input$simulate == 0 || is.null(values[["sim"]])){return()}
      
      if(input$pre_elast == "Pre-Merger"){ preMerger = TRUE}
      else{preMerger =FALSE}
      
      res <- as.matrix(elast(values[["sim"]], preMerger=preMerger, market = TRUE))
      colnames(res)= "Market"
      print(res)
      
    }, rownames = FALSE)   


    ## display messages to message tab
   output$warnings <- renderPrint({
    
      if(input$inTabset!= "msgpanel" || input$simulate == 0 || is.null(values[["msg"]]$warning)){return()}  
     
      print(values[["msg"]]$warning)
      
      
   })  

   output$errors <- renderPrint({
     
     if(input$inTabset!= "msgpanel" || input$simulate == 0 || is.null(values[["msg"]]$error)){cat(return())}  
     
     print(values[["msg"]]$error)
     
     
   })  

  })

