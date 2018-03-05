require(shiny)
require(rhandsontable)


shinyServer(function(input, output, session) {


   nPossProds <- 10 #only allow 10 products
   
   
   genInputData <- function(){
    # a function to generate default input data set for simulations
     
     inputData <- data.frame(
       Name = c("Prod1","Prod2","Prod3","Prod4"),
       ownerPre  = c("Firm1","Firm2","Firm3","Firm3"),
       ownerPost = c("Firm1","Firm1","Firm3","Firm3"),
       'Prices \n(\u00A4/unit)'    =c(.0441,.0328,.0409,.0396)*100,
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
     
     thiscmcr <- NA
     try(thiscmcr <- cmcr(res), silent=TRUE)
     thiscv <- CV(res)
     thiselast <- elast(res,market=TRUE)
     
     
     thispsdelta  <- sum(calcProducerSurplus(res,preMerger=FALSE) -calcProducerSurplus(res,preMerger=TRUE),na.rm=TRUE)
     
     res <- data.frame(
           'HHI Change' = hhi(res,preMerger=FALSE) -  hhi(res,preMerger=TRUE),
           'Industry Price Change (%)' = sum((s$priceDelta)*s$sharesPost/100,na.rm=TRUE),
           'Merging Party Price Change (%)'= sum(((s$priceDelta)*s$sharesPost)[s$isParty == "*"] / s$sharesPost[s$isParty== "*"],na.rm=TRUE),
           'Compensating Marginal Cost Reduction (%)' = ifelse(isCournot, thiscmcr, sum(thiscmcr*s$sharesPost[s$isParty== "*"]/sum(s$sharesPost[s$isParty== "*"],na.rm=TRUE))),
           'Consumer Harm (\u00A4/unit)' = -1*thiscv,
           'Producer Benefit (\u00A4/unit)' = thispsdelta,
           'Overall Effect (\u00A4/unit)'= -1*thiscv + thispsdelta,
           
           'Estimated Market Elasticity' = thiselast,
           check.names=FALSE
     )
     
     if(isRevDemand){
       res$'Overall Effect (\u00A4/unit)' <- NULL
       colnames(res) <- gsub('(?<=Consumer Harm\\s)\\(\\\u00A4/unit\\)',"(% Expenditure)",colnames(res), perl=TRUE)
       }
     
     
     if(is.na(res[,"Compensating Marginal Cost Reduction (%)"])) res[,"Compensating Marginal Cost Reduction (%)"] <- NULL
     
     return(res)
   }


   
   runSims <- function(supply, demand, indata, mktElast){
     # a function to execute code from antitrust package based on ui inputs
     
     firstMargin <- which(!is.na(indata$Margins))[1]
     firstPrice <- which(!is.na(indata[,"Prices \n(\u00A4/unit)"]))[1]
     
     
     ownerPre = model.matrix(~-1+indata$ownerPre)
     ownerPre = tcrossprod(ownerPre)
     ownerPost = model.matrix(~-1+indata$ownerPost)
     ownerPost = tcrossprod(ownerPost)
     
     
     switch(supply,
            Bertrand =
              switch(demand,
                     `logit (unknown elasticity)`= logit.alm(prices= indata[,"Prices \n(\u00A4/unit)"],
                                                             shares= indata[,"Output"],
                                                             margins= indata$Margins,
                                                             ownerPre= ownerPre,
                                                             ownerPost= ownerPost,
                                                             mcDelta = indata$mcDelta, labels=indata$Name),
                     `ces (unknown elasticity)`= ces.alm(prices= indata[,"Prices \n(\u00A4/unit)"],
                                                         shares= indata[,"Output"],
                                                         margins= indata$Margins,
                                                         ownerPre= ownerPre,
                                                         ownerPost= ownerPost,
                                                         mcDelta = indata$mcDelta, labels=indata$Name),
                     linear=linear(prices= indata[,"Prices \n(\u00A4/unit)"],
                                   quantities= indata[,"Output"],
                                   margins= indata$Margins,
                                   ownerPre= ownerPre,
                                   ownerPost= ownerPost,
                                   mcDelta = indata$mcDelta, labels=indata$Name),
                     aids=aids(prices= indata[,"Prices \n(\u00A4/unit)"],
                               shares= indata[,"Output"],
                               margins= indata$Margins,
                               ownerPre= ownerPre,
                               ownerPost= ownerPost,
                               mcDelta = indata$mcDelta, labels=indata$Name),
                     logit= logit.alm(prices= indata[,"Prices \n(\u00A4/unit)"],
                                      shares= indata[,"Output"],
                                      margins= indata$Margins,
                                      ownerPre= ownerPre,
                                      ownerPost= ownerPost,
                                      mcDelta = indata$mcDelta, labels=indata$Name,  mktElast = mktElast ),
                     ces = ces.alm(prices= indata[,"Prices \n(\u00A4/unit)"],
                                   shares= indata[,"Output"],
                                   margins= indata$Margins,
                                   ownerPre= ownerPre,
                                   ownerPost= ownerPost,
                                   mcDelta = indata$mcDelta, labels=indata$Name,  mktElast = mktElast),
                     linear=linear(prices= indata[,"Prices \n(\u00A4/unit)"],
                                   quantities= indata[,"Output"],
                                   margins= indata$Margins,
                                   ownerPre= ownerPre,
                                   ownerPost= ownerPost,
                                   mcDelta = indata$mcDelta, labels=indata$Name),
                     pcaids=pcaids(prices= indata[,"Prices \n(\u00A4/unit)"],
                                   shares= indata[,"Output"],
                                   knownElast = -1/indata$Margins[firstMargin],
                                   knownElastIndex = firstMargin,
                                   mktElast = mktElast,
                                   ownerPre= ownerPre,
                                   ownerPost= ownerPost,
                                   mcDelta = indata$mcDelta, labels=indata$Name)
              ),
            Cournot = 
              
              cournot(prices= indata[,"Prices \n(\u00A4/unit)"][firstPrice],
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
                                        `logit (unknown elasticity)` = auction2nd.logit.alm(prices= indata[,"Prices \n(\u00A4/unit)"],
                                                                                            shares= indata[,"Output"],
                                                                                            margins= indata$Margins,
                                                                                            ownerPre= ownerPre,
                                                                                            ownerPost= ownerPost,
                                                                                            mcDelta = indata$mcDelta, labels=indata$Name),
                                        logit = auction2nd.logit.alm(prices= indata[,"Prices \n(\u00A4/unit)"],
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
   values <- reactiveValues(inputData = genInputData(), sim = NULL)
      
   
                        
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
      

     updateTabsetPanel(session,inputId  = "inTabset", selected = "respanel")
      
      indata <- values[["inputData"]]
      
      colnames(indata)[colnames(indata) %in% c("Quantities","Shares")] <- "Output"
      
      indata <- indata[!is.na(indata[,"Output"]),]
      
      
      
      indata$mcDelta <- 0
      
      indata$ownerPre <- factor(indata$ownerPre)
      indata$ownerPost <- factor(indata$ownerPost,levels=levels(indata$ownerPre))
      
      
     
      

     values[["sim"]] <-  runSims(supply = input$supply,demand = demand(), indata = indata, mktElast = input$enterElast )
               
               
      
      
    })
    
    
   
      
    ## display summary results from gensum to results tab  
    output$results <-
          
          renderTable({
           
            if(input$inTabset != "respanel" || input$simulate == 0){return()}
         
           
            gensum(values[["sim"]])
          })
  
    
  ## display summary values to details tab      
    output$results_detailed <- renderTable({
       
      if(input$inTabset!= "detpanel" || input$simulate == 0 ){return()}
        
        if(input$supply == "Cournot"){
          
          
          capture.output(res <- summary(values[["sim"]], market=FALSE))
          
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
      
      if(input$inTabset!= "elastpanel" || input$simulate == 0 ){return()}
        
       if(input$pre_elast == "Pre-Merger"){ preMerger = TRUE}
        else{preMerger =FALSE}
        
        res <- elast(values[["sim"]], preMerger=preMerger)
        print(res)
        
      
    }, rownames = TRUE)   
  





  })

