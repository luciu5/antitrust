require(shiny)
require(rhandsontable)


shinyServer(function(input, output, session) {


   nPossProds <- 10 #only allow 10 products
   
   
   genInputData <- function(){

     
     inputData <- data.frame(
       Name = c("P1","P2","P3","P4"),
       ownerPre  = c("O1","O2","O3","O3"),
       ownerPost = c("O1","O1","O3","O3"),
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
           'Market Elasticity' = thiselast,
           check.names=FALSE
     )
     
     if(demand %in% c("aids","ces")){
       res$'Overall Effect ($/unit)' <- NULL
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
      demand <- ifelse(supply =="Bertrand", input$demand_bert,
                       ifelse(supply == "Cournot", input$demand_cournot, input$demand_2nd))

      updateCheckboxInput(session, inputId = "dispDetails", value = FALSE)      

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
      demand <- ifelse(supply =="Bertrand", input$demand_bert,
                       ifelse(supply == "Cournot", input$demand_cournot, input$demand_2nd))
       
     
      #if(input$dispDetails){updateCheckboxInput(session, inputId = "dispDetails", value = TRUE)}
      
      indata <- values[["inputData"]]
      
      outstring <- colnames(indata[colnames(indata) %in% c("Quantities","Shares")])
      
      indata <- indata[!is.na(indata[,outstring]),]
      
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
                        logit= logit.alm(prices= indata[,"Prices (\u00A4)"],
                                      shares= indata[,outstring],
                                      margins= indata$Margins,
                                      ownerPre= ownerPre,
                                      ownerPost= ownerPost,
                                      mcDelta = indata$mcDelta, labels=indata$Name),
                        ces= ces.alm(prices= indata[,"Prices (\u00A4)"],
                                 shares= indata[,outstring],
                                 margins= indata$Margins,
                                 ownerPre= ownerPre,
                                 ownerPost= ownerPost,
                                 mcDelta = indata$mcDelta, labels=indata$Name),
                        # linear=linear(prices= indata[,"Prices (\u00A4)"],
                        #               quantities= indata[,outstring],
                        #               margins= indata$Margins,
                        #               ownerPre= ownerPre,
                        #               ownerPost= ownerPost,
                        #               mcDelta = indata$mcDelta, labels=indata$Name),
                        aids=aids(prices= indata[,"Prices (\u00A4)"],
                                  shares= indata[,outstring],
                                  margins= indata$Margins,
                                  ownerPre= ownerPre,
                                  ownerPost= ownerPost,
                                  mcDelta = indata$mcDelta, labels=indata$Name)
                 ),
               Cournot =                           cournot(prices= indata[,"Prices (\u00A4)"][1],
                                                           demand = demand,
                                                           cost= rep("linear", nrow(indata)),
                                                           quantities = as.matrix(indata[,outstring]),
                                                           margins= as.matrix(indata$Margins),
                                                           ownerPre= ownerPre,
                                                           ownerPost= ownerPost,
                                                           mcDelta = indata$mcDelta, 
                                                           labels=list(as.character(indata$ownerPre),indata$Name[1])),
               `2nd Score Auction`=auction2nd.logit.alm(prices= indata[,"Prices (\u00A4)"],
                                                    shares= indata[,outstring],
                                                    margins= indata$Margins,
                                                    ownerPre= ownerPre,
                                                    ownerPost= ownerPost,
                                                    mcDelta = indata$mcDelta, labels=indata$Name)
               
               
        )
      
      
    })
    })
    
    
    thisSim <- NULL
    
    output$results <-
          
          renderTable({
           
            if(is.null(input$hot)){return()}
            
            thisSim <<- sims()
            supply <- input$supply
            demand <- ifelse(supply =="Bertrand", input$demand_bert,
                             ifelse(supply == "Cournot", input$demand_cournot, input$demand_2nd))
            
            #if(input$simulate == 0) return(cat("Enter Inputs"))
            gensum(thisSim,supply,demand)
          })
        
    output$results_detailed <- renderTable({
       
      if(input$dispDetails){
        
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
    
    
    
  





  })

