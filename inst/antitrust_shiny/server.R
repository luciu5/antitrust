require(shiny)
require(rhandsontable)


shinyServer(function(input, output, session) {


   nPossProds <- 10 #only allow 10 products
   
   
   genInputData <- function(){

     
     inputData <- data.frame(
       Name = c("Prod1","Prod2","Prod3","Prod4"),
       ownerPre  = c("Firm1","Firm2","Firm3","Firm3"),
       ownerPost = c("Firm1","Firm1","Firm3","Firm3"),
       #'Prices (\u00A4)'  = c(6.46, 6.25, 6.10),
       #Shares = c( 0.2643934, 0.3438824, 0.3917243),
       #Margins = c(0.3, 0.4,0.28),
       'Prices (\u00A4)'    =c(.0441,.0328,.0409,.0396)*100,
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
   
   gensum <- function(res,supply,demand){
     
     if(supply=="Cournot"){ 
      
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
           'Consumer Harm (\u00A4/unit)' = -1*thiscv,
           'Producer Benefit (\u00A4/unit)' = thispsdelta,
           'Overall Effect (\u00A4/unit)'= -1*thiscv + thispsdelta,
           'Compensating Marginal Cost Reduction (\u00A4/unit)' = ifelse(supply == "Cournot", thiscmcr, sum(thiscmcr*s$sharesPost[s$isParty== "*"]/sum(s$sharesPost[s$isParty== "*"],na.rm=TRUE))),
           'Estimated Market Elasticity' = thiselast,
           check.names=FALSE
     )
     
     if(demand %in% c("aids","ces","pcaids")){
       res$'Overall Effect (\u00A4/unit)' <- NULL
       colnames(res) <- gsub('(?<=Consumer Harm\\s)\\(\\\u00A4/unit\\)',"(% Expenditure)",colnames(res), perl=TRUE)
       }
     
     
     if(is.na(res[,"Compensating Marginal Cost Reduction (\u00A4/unit)"])) res[,"Compensating Marginal Cost Reduction (\u00A4/unit)"] <- NULL
     
     return(res)
   }


   values <- reactiveValues(inputData = NULL)


    observe({

      
     values$inputData <- genInputData()

    })

    observe({

       supply <- input$supply
       calcElast <- input$calcElast
       mktElast <- input$enterElast
       demand <- ifelse(supply =="Bertrand", 
                        ifelse(calcElast,input$demand_bert_alm,input$demand_bert),
                        ifelse(supply == "Cournot", 
                               ifelse(calcElast,input$demand_cournot_alm, input$demand_cournot), 
                               ifelse(calcElast,input$demand_2nd_alm, input$demand_2nd)))
       
       
      
      #updateCheckboxInput(session, inputId = "dispDetails", value = FALSE)      
      updateTabsetPanel(session,inputId  = "inTabset", selected = "respanel")
      
      if (!is.null(input$hot) 
          #&& input$simulate == 0
          ) {
        values$inputData = hot_to_r(input$hot)
      }
    })
    

    output$hot <- renderRHandsontable({
      
      inputData <- values[["inputData"]]
      
      if(input$supply =="Cournot"){colnames(inputData)[colnames(inputData) == "Shares"] <- "Quantities"}
      else{colnames(inputData)[colnames(inputData) == "Quantities"] <- "Shares"}
      
      if (!is.null(inputData))
        rhandsontable(inputData, stretchH = "all")
    })


    
    sims <- reactive({
      
      
      if(is.null(input$hot)){return()}
      
      supply <- input$supply
      calcElast <- input$calcElast
      demand <- ifelse(supply =="Bertrand", 
                       ifelse(calcElast,input$demand_bert_alm,input$demand_bert),
                       ifelse(supply == "Cournot", 
                              ifelse(calcElast,input$demand_cournot_alm, input$demand_cournot), 
                              ifelse(calcElast,input$demand_2nd_alm, input$demand_2nd)))
      
      mktElast <- input$enterElast
      
      
      
      #if(input$dispDetails){updateCheckboxInput(session, inputId = "dispDetails", value = TRUE)}
      
      indata <- values[["inputData"]]
      
      outstring <- colnames(indata[colnames(indata) %in% c("Quantities","Shares")])
      
      indata <- indata[!is.na(indata[,outstring]),]
      
      
      firstMargin <- which(!is.na(indata$Margins))[1]
      firstPrice <- which(!is.na(indata[,"Prices (\u00A4)"]))[1]
      
      #if(!input$incEff) indata$mcDelta <- 0
      indata$mcDelta <- 0
      
      indata$ownerPre <- factor(indata$ownerPre)
      indata$ownerPost <- factor(indata$ownerPost,levels=levels(indata$ownerPre))
      
      
      ownerPre = model.matrix(~-1+indata$ownerPre)
      ownerPre = tcrossprod(ownerPre)
      ownerPost = model.matrix(~-1+indata$ownerPost)
      ownerPost = tcrossprod(ownerPost)
    
    isolate({  
     
        switch(supply,
               Bertrand =
                 switch(demand,
                        `logit (unknown elasticity)`= logit.alm(prices= indata[,"Prices (\u00A4)"],
                                      shares= indata[,outstring],
                                      margins= indata$Margins,
                                      ownerPre= ownerPre,
                                      ownerPost= ownerPost,
                                      mcDelta = indata$mcDelta, labels=indata$Name),
                        `ces (unknown elasticity)`= ces.alm(prices= indata[,"Prices (\u00A4)"],
                                 shares= indata[,outstring],
                                 margins= indata$Margins,
                                 ownerPre= ownerPre,
                                 ownerPost= ownerPost,
                                 mcDelta = indata$mcDelta, labels=indata$Name),
                         linear=linear(prices= indata[,"Prices (\u00A4)"],
                                       quantities= indata[,outstring],
                                       margins= indata$Margins,
                                       ownerPre= ownerPre,
                                       ownerPost= ownerPost,
                                       mcDelta = indata$mcDelta, labels=indata$Name),
                        aids=aids(prices= indata[,"Prices (\u00A4)"],
                                  shares= indata[,outstring],
                                  margins= indata$Margins,
                                  ownerPre= ownerPre,
                                  ownerPost= ownerPost,
                                  mcDelta = indata$mcDelta, labels=indata$Name),
                        logit= logit.alm(prices= indata[,"Prices (\u00A4)"],
                                     shares= indata[,outstring],
                                     margins= indata$Margins,
                                     ownerPre= ownerPre,
                                     ownerPost= ownerPost,
                                     mcDelta = indata$mcDelta, labels=indata$Name,  mktElast = mktElast ),
                        ces = ces.alm(prices= indata[,"Prices (\u00A4)"],
                                 shares= indata[,outstring],
                                 margins= indata$Margins,
                                 ownerPre= ownerPre,
                                 ownerPost= ownerPost,
                                 mcDelta = indata$mcDelta, labels=indata$Name,  mktElast = mktElast),
                        linear=linear(prices= indata[,"Prices (\u00A4)"],
                                      quantities= indata[,outstring],
                                      margins= indata$Margins,
                                      ownerPre= ownerPre,
                                      ownerPost= ownerPost,
                                      mcDelta = indata$mcDelta, labels=indata$Name),
                        pcaids=pcaids(prices= indata[,"Prices (\u00A4)"],
                                      shares= indata[,outstring],
                                      knownElast = -1/indata$Margins[firstMargin],
                                      knownElastIndex = firstMargin,
                                      mktElast = mktElast,
                                      ownerPre= ownerPre,
                                      ownerPost= ownerPost,
                                      mcDelta = indata$mcDelta, labels=indata$Name)
                 ),
               Cournot = 

                                 cournot(prices= indata[,"Prices (\u00A4)"][firstPrice],
                                                                                   demand = gsub("\\s+\\(.*","",demand,perl=TRUE),
                                                                                   cost= rep("linear", nrow(indata)),
                                                                                   quantities = as.matrix(indata[,outstring]),
                                                                                   margins= as.matrix(indata$Margins),
                                                                                   ownerPre= ownerPre,
                                                                                   ownerPost= ownerPost,
                                                                                   mktElast = ifelse( grepl("unknown elasticity", demand),
                                                                                                      NA_real_, mktElast),
                                                                                   mcDelta = indata$mcDelta, 
                                                                                   labels=list(as.character(indata$ownerPre),indata$Name[firstPrice]))
                                ,
               `2nd Score Auction`= switch(demand,
                                           `logit (unknown elasticity)` = auction2nd.logit.alm(prices= indata[,"Prices (\u00A4)"],
                                                    shares= indata[,outstring],
                                                    margins= indata$Margins,
                                                    ownerPre= ownerPre,
                                                    ownerPost= ownerPost,
                                                    mcDelta = indata$mcDelta, labels=indata$Name),
                                            logit = auction2nd.logit.alm(prices= indata[,"Prices (\u00A4)"],
                                                                                shares= indata[,outstring],
                                                                                margins= indata$Margins,
                                                                                ownerPre= ownerPre,
                                                                                ownerPost= ownerPost,
                                                                                mcDelta = indata$mcDelta, labels=indata$Name, 
                                                                               mktElast = mktElast)
                                           
               )
               
               
        )
      
      
    })
        
      
      
     
      
      
    })
    
    
    thisSim <- NULL
    
    output$results <-
          
          renderTable({
           
            if(is.null(input$hot)){return()}
         
            if(input$inTabset!= "respanel"){return()}   
            thisSim <<- sims()
            supply <- input$supply
            calcElast <- input$calcElast
            demand <- ifelse(supply =="Bertrand", 
                             ifelse(calcElast,input$demand_bert_alm,input$demand_bert),
                             ifelse(supply == "Cournot", 
                                    ifelse(calcElast,input$demand_cournot_alm, input$demand_cournot), 
                                    ifelse(calcElast,input$demand_2nd_alm, input$demand_2nd)))
            
            
            
            #if(input$simulate == 0) return(cat("Enter Inputs"))
            gensum(thisSim,supply,demand)
          })
        
    output$results_detailed <- renderTable({
       
      if(input$inTabset== "detpanel"){
        
        if(input$supply == "Cournot"){
          
          
          capture.output(res <- summary(thisSim, market=FALSE))
          
          res$product <- res$mcDelta <- NULL
          
          try(colnames(res) <- c("Party","Name","Pre-Merger Price","Post-Merger Price", "Price Change (%)","Pre-Merger Quantity","Post-Merger Quantity", "Output Change (%)"),silent=TRUE)
          
          }
        
        else{
          
          capture.output(res <- summary(thisSim))
          res$Name <- rownames(res)
          res <- res[,c(1, ncol(res), 2 : (ncol(res)  - 1))]
          colnames(res) <- c("Party","Name","Pre-Merger Price","Post-Merger Price", "Price Change (%)","Pre-Merger Share (%)","Post-Merger Share (%)", "Share Change (%)")
        }
        
        
        print(res)
        
      }
      })  
    
    
   
    
    output$results_elast <- renderTable({
      
      if(input$inTabset== "elastpanel"){
        
       if(input$pre_elast == "Pre-Merger"){ preMerger = TRUE}
        else{preMerger =FALSE}
        
        res <- elast(thisSim, preMerger=preMerger)
        print(res)
        
      }
    }, rownames = TRUE)   
  





  })

