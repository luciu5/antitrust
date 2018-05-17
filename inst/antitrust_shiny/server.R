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
       'Pre-merger\n Owner'  = c("Firm1","Firm2","Firm3","Firm3"),
       'Post-merger\n Owner' = c("Firm1","Firm1","Firm3","Firm3"),
       'Prices \n($/unit)'    = rep(10,4),
       'Quantities'   =c(0.4,.3,.2,.1)*100,
       'Margins\n(p-c)/p' =c(0.5,NA,NA,NA),
       'Post-merger\n Cost Changes\n(Proportion)' = rep(0,4),
       stringsAsFactors = FALSE,
       check.names=FALSE
     )
     
     
     
     nDefProd <- nrow(inputData)
     inputData <- inputData[c(1:nDefProd,rep(1, nPossProds - nDefProd)),]
     #if(input$incEff) inputData$mcDelta <- 0
     
     inputData[(nDefProd+1):nPossProds,] <- NA
     inputData <- inputData[order(inputData$`Quantities`, decreasing = TRUE),]
     rownames(inputData) <- NULL
     
    
     return(inputData)

   }
   
   gensum <- function(res){
     #a function to generate summary stats for results tab
     
       
     isCournot <- grepl("Cournot",class(res))
     isAuction <- grepl("Auction",class(res))
     isRevDemand <- grepl("ces|aids",class(res),ignore.case = TRUE)
     isLogit <- grepl("logit",class(res),ignore.case = TRUE)
     
     missPrices <- any(is.na(res@prices))
     
     inLevels <- FALSE
     
     if(isAuction && missPrices){inLevels = TRUE}

     
     
     
     if(isCournot){ 
      
       capture.output( s <- summary(res, market = FALSE))
       theseshares <- drop(res@quantities/sum(res@quantities))
       
       totQuantPost <- sum(s$quantityPost,na.rm=TRUE)
       s$sharesPost <- s$quantityPost/totQuantPost*100
       
     }
     
     else{
       capture.output(s <- summary(res, revenue = isRevDemand & missPrices,levels = inLevels))
       theseshares <- calcShares(res, preMerger=TRUE, revenue=isRevDemand & missPrices)
       
       
       theseshares <- theseshares/sum(theseshares)
       
       }
     
     isparty <- s$isParty == "*"
     
     thispsdelta <- thiscmcr <- thiscv <- NA
     
     try(thiscmcr <- cmcr(res), silent=TRUE)
     try(thiscv <- CV(res),silent = TRUE)
     
     
     try(thispsdelta  <- sum(calcProducerSurplus(res,preMerger=FALSE) - calcProducerSurplus(res,preMerger=TRUE)),silent=TRUE)
     
     
     partyshare <- s$sharesPost[isparty]
     
     
     result <- with(s, data.frame(
           #'HHI Change' = as.integer(round(hhi(res,preMerger=FALSE) -  hhi(res,preMerger=TRUE))),
           'Pre-Merger HHI' =  as.integer(HHI(theseshares,owner=res@ownerPre)),
           'HHI Change' = as.integer(HHI(theseshares,owner=res@ownerPost) - HHI(theseshares,owner=res@ownerPre)),
           'Industry Price Change (%)' = sum(priceDelta * sharesPost/100,na.rm=TRUE),
           'Merging Party Price Change (%)'= sum(priceDelta[isparty] * partyshare, na.rm=TRUE) / sum(partyshare),
           'Compensating Marginal Cost Reduction (%)' = ifelse(isCournot, thiscmcr, sum(thiscmcr * partyshare) / sum(partyshare)),
           'Consumer Harm ($)' = thiscv,
           'Producer Benefit ($)' = thispsdelta,
           'Difference ($)'= thiscv - thispsdelta,
           
           #'Estimated Market Elasticity' = thiselast,
           check.names=FALSE
     ))
    
     
     ## Append parenthetical with % of post-merger revenues
     result <- rbind(result,result)
     result[2, !grepl("\\$",colnames(result))] <- NA
     result[,-(1:2)] <- round(result[,-(1:2)],digits = 1)
     result[2,grepl("\\$",colnames(result))] <- paste0("(",round(result[2 , grepl("\\$",colnames(result))]*100 / calcRevenues(res, preMerger = FALSE, market=TRUE), digits = 1), "%)")
     
     
     if(inLevels){
       colnames(result) <- gsub('(?<=Price Change\\s)\\(%\\)',"($/unit)",colnames(result), perl=TRUE)
     }
     
    
    
 
     if(all(is.na(result[,"Compensating Marginal Cost Reduction (%)"]))) result[,"Compensating Marginal Cost Reduction (%)"] <- NULL
     
     return(result)
   }


   
   gendiag <- function(res,mktElast=FALSE){
     #a function to generate diagnostics data
     
    isCournot <- grepl("Cournot",class(res))
    
    if(isCournot){labels= res@labels[[1]]}
    else{labels=res@labels}
    
    obsPrices <- res@prices
    obsShares <- res@shares
    obsMargins <- res@margins
    obsElast <- res@mktElast
    
    #if(length(obsMargins[!is.na(obsMargins)]) < 2){return()}
    
    
    prePrices <- unname(drop(res@pricePre))
    preMargins <- drop(calcMargins(res, preMerger=TRUE))
    preShares <- drop(calcShares(res, preMerger=TRUE))
    preShares <- drop(preShares/sum(preShares))
    preElast <- elast(res, preMerger=TRUE, market=TRUE)
    
    if(!mktElast){
      
    res <- data.frame(
      "Inputted Prices"= obsPrices,
      "Fitted Prices" = prePrices,
      "Price Change (%)"= (1 - obsPrices/prePrices)*100,
      "Inputted Shares (%)" = obsShares*100,
      "Fitted Shares(%)"=preShares*100,
      "Share Change (%)"=(1 - obsShares/preShares)*100,
      "Inputted Margins" = obsMargins,
      "Fitted  Margins"=preMargins,
      "Margin Change (%)"= (1 - obsMargins/preMargins)*100,
      #'Market Elasticity'= 1 - obsElast/preElast,
      check.names = FALSE
    )
    
    #rmThese <- colSums(abs(res),na.rm=TRUE)
    
  
   
    if(isCournot)  res[-1,grepl('Prices',colnames(res))] <- NA
     
    #res <- res[,rmThese >1e-3,drop=FALSE]
    
    
    if(!isCournot) rownames(res) <- labels
    
    }
    
    else{ res <- data.frame(
                    'Inputted Elasticity' = obsElast,
                    'Fitted Elasticity' = preElast,
                     'Elasticity Change'= (1 - obsElast/preElast)*100,
                            check.names = FALSE)
    
    #if(res < 1e-3) res <- NULL
    }
    return(res)
   }
   
   runSims <- function(supply, demand, indata, mktElast){
     # a function to execute code from antitrust package based on ui inputs
     
     prices <- indata[,"Prices \n($/unit)"]
     margins <- indata$Margins
     
     missPrices <- any(is.na(prices))
     
     shares_quantity <- shares_revenue <- indata$Output/sum(indata$Output, na.rm=TRUE)
     insideSize <- sum(indata[,"Output"], na.rm=TRUE)
     
     
     if(!missPrices){ 
       
       if(any(grepl("ces|aids",demand,perl=TRUE,ignore.case=TRUE))) insideSize <- sum(prices * indata[,"Output"], na.rm =TRUE)
       
       shares_revenue <- prices * shares_revenue / sum(prices * shares_revenue)
       if(supply == "2nd Score Auction") margins <- margins * prices # convert to level margins
     }
     
     firstMargin <- which(!is.na(margins))[1]
     firstPrice <- which(!is.na(prices))[1]
     
     
     ownerPre = model.matrix(~-1+indata$'Pre-merger\n Owner')
     ownerPre = tcrossprod(ownerPre)
     ownerPost = model.matrix(~-1+indata$'Post-merger\n Owner')
     ownerPost = tcrossprod(ownerPost)
     
     
     
     
     switch(supply,
            Bertrand =
              switch(demand,
                     `logit (unknown elasticity)`= logit.alm(prices= prices,
                                                             shares= shares_quantity,
                                                             margins= margins,
                                                             ownerPre= ownerPre,
                                                             ownerPost= ownerPost,
                                                             insideSize = insideSize ,
                                                             mcDelta = indata$mcDelta, labels=indata$Name),
                     `aids (unknown elasticity)` = aids(prices= prices,
                               shares= shares_revenue,
                               margins= margins,
                               ownerPre= ownerPre,
                               ownerPost= ownerPost,
                               insideSize = insideSize ,
                               mcDelta = indata$mcDelta, labels=indata$Name),
                     `ces (unknown elasticity)`= ces.alm(prices= prices,
                                                         shares= shares_revenue,
                                                         margins= margins,
                                                         ownerPre= ownerPre,
                                                         ownerPost= ownerPost,
                                                         insideSize = insideSize ,
                                                         mcDelta = indata$mcDelta, labels=indata$Name),
                     linear=linear(prices= prices,
                                   quantities= indata[,"Output"],
                                   margins= margins,
                                   ownerPre= ownerPre,
                                   ownerPost= ownerPost,
                                   mcDelta = indata$mcDelta, labels=indata$Name),
                     aids=aids(prices= prices,
                               shares= shares_revenue,
                               margins= margins,
                               ownerPre= ownerPre,
                               ownerPost= ownerPost,
                               insideSize = insideSize ,
                               mcDelta = indata$mcDelta, labels=indata$Name, mktElast = mktElast),
                     logit= logit.alm(prices= prices,
                                      shares= shares_quantity,
                                      margins= margins,
                                      ownerPre= ownerPre,
                                      ownerPost= ownerPost,
                                      insideSize = insideSize ,
                                      mcDelta = indata$mcDelta, labels=indata$Name,  mktElast = mktElast ),
                     ces = ces.alm(prices= prices,
                                   shares= shares_revenue,
                                   margins= margins,
                                   ownerPre= ownerPre,
                                   ownerPost= ownerPost,
                                   insideSize = insideSize ,
                                   mcDelta = indata$mcDelta, labels=indata$Name,  mktElast = mktElast),
                     linear=linear(prices= prices,
                                   quantities= indata[,"Output"],
                                   margins= margins,
                                   ownerPre= ownerPre,
                                   ownerPost= ownerPost,
                                   mcDelta = indata$mcDelta, labels=indata$Name),
                     pcaids=pcaids(prices= prices,
                                   shares= shares_revenue,
                                   knownElast = -1/margins[firstMargin],
                                   knownElastIndex = firstMargin,
                                   mktElast = mktElast,
                                   insideSize = insideSize,
                                   ownerPre= ownerPre,
                                   ownerPost= ownerPost,
                                   mcDelta = indata$mcDelta, labels=indata$Name)
              ),
            Cournot = 
              
              cournot(prices= prices[firstPrice],
                      demand = gsub("\\s+\\(.*","",demand,perl=TRUE),
                      cost= rep("linear", nrow(indata)),
                      quantities = as.matrix(indata[,"Output"]),
                      margins= as.matrix(margins),
                      ownerPre= ownerPre,
                      ownerPost= ownerPost,
                      mktElast = ifelse( grepl("unknown elasticity", demand),
                                         NA_real_, mktElast),
                      mcDelta = indata$mcDelta, 
                      labels=list(as.character(indata$Name),indata$Name[firstPrice]))
            ,
            `2nd Score Auction`= switch(demand,
                                        `logit (unknown elasticity)` = auction2nd.logit.alm(prices= prices,
                                                                                            shares= shares_quantity,
                                                                                            margins= margins,
                                                                                            ownerPre= ownerPre,
                                                                                            ownerPost= ownerPost,
                                                                                            mcDelta = indata$mcDelta, labels=indata$Name),
                                        logit = auction2nd.logit.alm(prices= prices,
                                                                     shares= shares_quantity,
                                                                     margins= margins,
                                                                     ownerPre= ownerPre,
                                                                     ownerPost= ownerPost,
                                                                     insideSize = insideSize,
                                                                     mcDelta = indata$mcDelta, labels=indata$Name, 
                                                                     mktElast = mktElast)
                                        
            )
     )
   }
   
   
   
   # ## create a reactive object that tracks which demand is being used
   # demand <- eventReactive(input$simulate, {
   #                                         ifelse(input$supply =="Bertrand", 
   # 
   #                                            ifelse(input$calcElast == "2 or more margins",input$demand_bert_alm,input$demand_bert),
   #                                                    ifelse(input$supply == "Cournot",
   #                                                   ifelse(input$calcElast == "1 or more margins",input$demand_cournot_alm, input$demand_cournot),
   #                                                   ifelse(input$calcElast == "2 or more margins",input$demand_2nd_alm, input$demand_2nd)))
   #                    })
   # 
   
   ## create a reactive list of objects
   values <- reactiveValues(inputData = genInputData(), sim = NULL, msg = NULL)
      
  
                        
   ## update reactive list whenever changes are made to input
    observe({
    
      supply <- input$supply
      demand <- gsub("\\s*(unknown elasticity)","",input$demand,perl=TRUE)
      
      provElast <- grepl('elasticity',input$calcElast)
      defElast <- ifelse(provElast, 1 , 2)
      
     if(supply == 'Cournot'){
       theseChoices <- c("market elasticity and 0 or more margins",
                         "1 or more margins"
       )
       
       
       demandChoices<- c("linear","log")
       
       }
      
    else{
      theseChoices <- c("market elasticity and 1 or more margins",
                        "2 or more margins"
      )
      updateRadioButtons(session = session, "calcElast", "Calibrate model parameters using:",
                         choices = theseChoices , selected = theseChoices[defElast]) 
      
      if(supply == 'Bertrand'){demandChoices<- c('logit','ces','aids')}
      else{demandChoices<- c('logit')}
      
     
      
    }
      
      if( !provElast ){ demandChoices <- paste(demandChoices,"(unknown elasticity)")}
      
      updateRadioButtons(session = session, "calcElast", "Calibrate model parameters using:",
                         choices = theseChoices , selected = theseChoices[defElast])
      
      provDemand <- grep(demand, demandChoices)
      provDemand <- ifelse(length(provDemand) > 0 , provDemand, 1)
      
      updateSelectInput(session=session, "demand", "Demand Specification:",
                        choices = demandChoices, selected = demandChoices[provDemand])
      
      
      if(!is.null(input$hot)){
        values$inputData = hot_to_r(input$hot)
        values$inputData[!is.na(values$inputData$Name) & values$inputData$Name != '',]
      }
    })
    

    ## display inputs 
    output$hot <- renderRHandsontable({
      
      inputData <- values[["inputData"]]
      
      prices <- inputData[,"Prices \n($/unit)"]
      output <- inputData[,grepl("Quantities|Revenue",colnames(inputData), perl=TRUE)]
      
      missPrices <- isTRUE(any(is.na(prices[ !is.na(output) ] ) ))
      
      if(input$supply == "2nd Score Auction"){colnames(inputData)[grepl("Cost Changes",colnames(inputData))] <-'Post-merger\n Cost Changes\n($/unit)'}
      else{colnames(inputData)[grepl("Cost Changes",colnames(inputData))] <-'Post-merger\n Cost Changes\n(Proportion)'}
      
      if(missPrices && input$supply =="2nd Score Auction"){colnames(inputData)[grepl("Margins",colnames(inputData))] <- "Margins\n ($/unit)"}
      else{colnames(inputData)[grepl("Margins",colnames(inputData))] <- "Margins\n (p-c)/p"}

      if (missPrices && any(grepl("ces|aids",input$demand, perl=TRUE), na.rm=TRUE)){colnames(inputData)[grepl("Quantities",colnames(inputData))] <- "Revenues"}
      else{{colnames(inputData)[grepl("Revenues",colnames(inputData))] <- "Quantities"}}
      
      if (!is.null(inputData))
        rhandsontable(inputData, stretchH = "all", contextMenu = FALSE ) %>% hot_col(col = 1:ncol(inputData), valign = "htMiddle") %>%
        hot_col(col = which (sapply(inputData,is.numeric)),halign = "htCenter" ) %>% hot_cols(columnSorting = TRUE)
    })


    
  ## simulate merger when the "simulate" button is clicked
   observeEvent(input$simulate,{
      

     values[["sim"]] <- values[["msg"]] <-  NULL
     
     updateTabsetPanel(session,inputId  = "inTabset", selected = "respanel")
      
      indata <- values[["inputData"]]
      
      isOutput <-  grepl("Quantities|Revenues",colnames(indata),perl = TRUE)
      
      
      colnames(indata)[ grepl("Margins",colnames(indata),perl = TRUE)] <- "Margins"
      
      colnames(indata)[isOutput] <- "Output"
      
      indata <- indata[!is.na(indata[,"Output"]),]
      
     
      
      indata$mcDelta <- indata[,grep('Cost Changes',colnames(indata))]
      indata$mcDelta[is.na(indata$mcDelta)] <- 0
      
      indata$'Pre-merger\n Owner' <- factor(indata$'Pre-merger\n Owner',levels=unique(indata$'Pre-merger\n Owner') )
      indata$'Post-merger\n Owner' <- factor(indata$'Post-merger\n Owner',levels=unique(indata$'Post-merger\n Owner'))
      
      
     
      thisSim <- msgCatcher(
                          runSims(supply = input$supply,demand = input$demand, indata = indata, mktElast = input$enterElast )
      )

    
     #thisSim$warning <- grep("Estimated outside share is close to 0", thisSim$warning, value= TRUE, invert=TRUE)
     #if(length(thisSim$warning) == 0){thisSim$warning = NULL}
       
     values[["sim"]] <-  thisSim$value
               
     values[["msg"]] <-  list(error=thisSim$error,warning=thisSim$warning)        
      
     if(!is.null(thisSim$error) || !is.null(thisSim$warning)) updateTabsetPanel(session,inputId  = "inTabset", selected = "msgpanel")
      
    })
    
  
   ## identify interface as Public Site or not
   output$overIDText <-   renderText({
     
     if(is.null(values[["inputData"]])){return()}
     
     provElast <- grepl('elasticity',input$calcElast)
     
     inputData <- values[["inputData"]]
     
     nMargins <- inputData[,grepl("Margins",colnames(inputData))]
     nMargins <- length(nMargins[!is.na(nMargins)])
     
     
     if(input$supply == "Cournot" && 
        ((provElast && nMargins > 0)  || (!provElast && nMargins >1) )){
       res <- paste(helpText(tags$b("Note:"), "some model parameters are over-identified. The tables above may be helpful in assessing model fit."))
     }
    
     else if(input$supply != "Cournot" && 
             ((provElast && nMargins > 1)  || (!provElast && nMargins >2) )){
       res <- paste(helpText(tags$b("Note:"),"some model parameters are over-identified. The tables above may be helpful in assessing model fit."))
     }
     else{
       res <- paste(helpText(tags$b("Note:"),"model parameters are just-identified. Inputted and fitted values should match."))
       
     }
     res
     
   })  
   
   ## identify interface as Public Site or not
   output$urlText <-   renderText({
     
     thisurl <- session$clientData$url_hostname
     
     
     if(grepl("atrnet\\.gov", thisurl, perl=TRUE)){
       res <- paste(h4(span("Internal Server",style="color:blue")))
     }
     else if(grepl("^127", thisurl, perl=TRUE)){
       res <- paste(h4(span("Local Server",style="color:green")))
     }
     else{
       res <- paste(h4(span("Public Server",style="color:red")))
     }
    res
    
   })
      
    ## display summary results from gensum to results tab  
    output$results <-
          
          renderTable({
           
            if(input$inTabset != "respanel" || input$simulate == 0|| is.null(values[["sim"]])){return()}
         
            inputData <- values[["inputData"]]
            
            gensum(values[["sim"]])
          }, na="", digits =1)
  
    
  ## display summary values to details tab      
    output$results_detailed <- renderTable({
       
      if(input$inTabset!= "detpanel" || input$simulate == 0  || is.null(values[["sim"]])){return()}
        
        if(input$supply == "Cournot"){
          
          res <- NULL
          capture.output(try(res <- summary(values[["sim"]], revenue= FALSE,market=FALSE),silent=TRUE))
          
          
          res$product <- res$mcDelta <- NULL
          
          try(colnames(res) <- c("Merging Party","Name","Pre-Merger Price","Post-Merger Price", "Price Change (%)","Pre-Merger Quantity","Post-Merger Quantity", "Output Change (%)"),silent=TRUE)
          
          }
        
        else{
          
          isAuction <- grepl("Auction",class(values[["sim"]]))
          isRevDemand <- grepl("ces|aids",class(values[["sim"]]),ignore.case = TRUE)
          inLevels <- FALSE
          #isAIDS <- grepl("aids",class(values[["sim"]]),ignore.case = TRUE)
          missPrice <- any(is.na(values[["sim"]]@prices)) 
          if(isAuction && missPrice){inLevels = TRUE}
          
          capture.output(res <- summary(values[["sim"]], revenue=isRevDemand & missPrice, insideOnly=TRUE, levels=inLevels))
          res$Name <- rownames(res)
          res$mcDelta <- NULL
          res <- res[,c(1, ncol(res), 2 : (ncol(res)  - 1))]
          res$cmcr <- NA
          try(res$cmcr[res$isParty=="*"] <- cmcr(values[["sim"]]))
          
          thesenames <-  c("Merging Party","Name","Pre-Merger Price","Post-Merger Price", "Price Change (%)","Pre-Merger Share (%)","Post-Merger Share (%)", "Share Change (%)",'Compensating Marginal Cost Reduction (%)')
          
          
          #if(isAIDS && missPrice){thesenames <- thesenames[!thesenames %in% c("Pre-Merger Price","Post-Merger Price")]}
          
          colnames(res) <- thesenames
          
          if(all(is.na(res$`Compensating Marginal Cost Reduction (%)`))) res$`Compensating Marginal Cost Reduction (%)` <- NULL
          
          #res[,c("Pre-Merger Share (%)","Post-Merger Share (%)")] <-  res[,c("Pre-Merger Share (%)","Post-Merger Share (%)")] * 100 / colSums( res[,c("Pre-Merger Share (%)","Post-Merger Share (%)")])
        
          if(inLevels){ colnames(res)[ colnames(res) == "Price Change (%)"] = "Price Change ($/unit)"}
          
          }
        
        
        res
        
      
      }, digits = 2)  
    
    
    output$results_shareOut <- renderTable({
      
      if(input$inTabset!= "detpanel" || input$simulate == 0  || is.null(values[["sim"]])){return()}
      if( grepl("cournot",class(values[["sim"]]),ignore.case = TRUE)){return()}
        
      isCES <- grepl("ces",class(values[["sim"]]),ignore.case = TRUE)
      
      res <- data.frame('No-purchase\n Share (%)'= c(
                         1 - sum(calcShares(values[["sim"]], preMerger=TRUE,revenue=isCES)),
                         1 - sum(calcShares(values[["sim"]], preMerger=FALSE,revenue=isCES))
                         )
                        ,check.names = FALSE
                        )*100
      
      res$'Revenues ($)' <- as.integer(round(c(calcRevenues(values[["sim"]], preMerger=TRUE, market = TRUE ),
                              calcRevenues(values[["sim"]], preMerger=FALSE, market = TRUE )
                              )))
    
      
      rownames(res) <- c("Pre-Merger","Post-Merger")
      
      if( grepl("aids",class(values[["sim"]]),ignore.case = TRUE)) res$'No-purchase\n Share (%)' <- NULL
      
      return(res)
      
      }, rownames = TRUE, digits=1,align="c")
    
    ## display results to diagnostics tab
    output$results_diagnostics <- renderTable({
      
      if(input$inTabset!= "diagpanel" || input$simulate == 0 || is.null(values[["sim"]])){return()}
      
      res <- gendiag(values[["sim"]])
      
      res
      
      
      
    }, digits = 0 ,rownames = TRUE,align="c")
    
    ## display market elasticity gap to diagnostics tab
    output$results_diag_elast <- renderTable({
      
      if(input$inTabset!= "diagpanel" || input$simulate == 0 || is.null(values[["sim"]])){return()}
      
      res <- gendiag(values[["sim"]], mktElast=TRUE)
      
      res
      
      
      
    }, digits =2,rownames = FALSE,align="c")
    
    
      
    ## display parameters to diagnostics tab
    output$parameters <- renderPrint({
      
      if(input$inTabset!= "diagpanel" || input$simulate == 0  || is.null(values[["sim"]])){return()}  
      
      print(getParms(values[["sim"]],digits=2))
      
      
    })  
    
   
    ## display elasticities to elasticity tab
    output$results_elast <- renderTable({
      
      if(input$inTabset!= "elastpanel" || input$simulate == 0 || is.null(values[["sim"]])){return()}
        
      isCournot <- grepl("Cournot",class(values[["sim"]]))
      
       if(input$pre_elast == "Pre-Merger"){ preMerger = TRUE}
        else{preMerger =FALSE}
        
      if(!isCournot && input$diversions){
        res <- diversion(values[["sim"]], preMerger=preMerger)
      }
      else{  res <- elast(values[["sim"]], preMerger=preMerger)}
        if(isCournot){colnames(res) <- "Elasticity"}
        
        res
        
      
    }, rownames = TRUE)   
  
    
    output$results_mktelast <- renderTable({
      
      if(input$inTabset!= "elastpanel" || input$simulate == 0 || is.null(values[["sim"]])){return()}
      
      if(input$pre_elast == "Pre-Merger"){ preMerger = TRUE}
      else{preMerger =FALSE}
      
      res <- as.matrix(elast(values[["sim"]], preMerger=preMerger, market = TRUE))
      colnames(res)= "Market"
      res
      
    }, rownames = FALSE)   


    ## display R code to code tab
    output$results_code <- renderPrint({
      
      if(input$inTabset!= "codepanel"){return()}  
      
       
      
      indata <- values[["inputData"]]
      indata <- indata[!is.na(indata$Name) & indata$Name != '',] 
      cat(NULL)
      
      cnames <- colnames(indata)
      cnames <- gsub("\n","",cnames)
      
      firstPrice <- which(!is.na(indata[,grep("price",cnames, ignore.case = TRUE)]))[1]
      
      argnames <- c("demand",
                    "ownerPre",
                    "ownerPost",
                    "prices",
                    "quantities",
                    "margins",
                    "mcDelta"
                    
                    )
      
      argvalues   <- paste0(argnames," = simdata$`",cnames,"`")
  
    thisElast <- ifelse(grepl("elast",input$calcElast), 
                        input$enterElast,
                        "NA_real_"
    )
    
    if(grepl("logit", input$demand, ignore.case =TRUE)){ 
                        thisSize <- "sum(simdata$`Quantities`)"
    }
    else if(all(is.na(indata[,grepl("price",cnames, ignore.case = TRUE)]))){
                              
      thisSize <-   "sum(simdata$`Revenues`)"
    }
    
    else{
      
      thisSize <- paste0("sum(simdata$`",grep("price",cnames, ignore.case = TRUE, value=TRUE),"`*simdata$`Quantities`)")
      
                       
    }
    
    thisdemand <- gsub("\\s*\\(.*","",input$demand,perl=TRUE)
      
      argvalues[1] <- paste(c("demand = ", shQuote(thisdemand)), collapse = "") 
     
      argvalues <- c(argvalues,
                     paste0("mktElast = ", thisElast,collapse = ""),
                     paste0("insideSize = ",thisSize, collapse=""),
                     "labels = simdata$Name"
                      
      )
      
      
      
      
      if(input$supply == "Cournot"){
        atrfun <- "cournot"
        argvalues[grep("prices", argvalues)] <- paste0(argvalues[grep("prices", argvalues)],"[",firstPrice,"]")
        argvalues[grep("quantities", argvalues)] <- "quantities = as.matrix(simdata$`Quantities`)" 
        argvalues[grep("margins", argvalues)] <- paste0("margins = as.matrix(simdata$`",grep("Margin",cnames,value = TRUE),"`)")

        argvalues[grep("labels", argvalues)] <- sprintf("labels = list(as.character(simdata$Name),as.character(simdata$Name[%d]))",firstPrice) 
        argvalues[grep("insideSize", argvalues)] <- NULL
        }
      else if( input$supply =="Bertrand"){atrfun <- "bertrand.alm"}
      else{atrfun <- "auction2nd.logit.alm"
           argvalues <- argvalues[-1]
           argvalues[grep("quantities", argvalues)] <- "shares = simdata$`Quantities` / sum( simdata$`Quantities` ) "
           argvalues[grep("margins", argvalues)] <- paste0("margins = simdata$`",
                                                           grep("Margin",cnames,value = TRUE),
                                                            "` * ", "simdata$`", grep("Price", cnames, value=TRUE),"`")
           }
      
      atrfun <- paste0("simres <- ",atrfun,"(\n\t",paste0(argvalues,collapse = ",\n\t"),")",collapse = "\n")
      
      indata_code <- sapply(1:ncol(indata),
                            function(x){d <- indata[,x];
                                      if(is.character(d)){d <- sprintf("'%s'", indata[,x])};
                                      paste0('`',cnames[x],'`',"= c(",paste(d,collapse=","),")")}
                            )
      indata_code <- paste0("simdata <- data.frame(\n\t", paste0(indata_code,collapse=",\n\t"),",\n check.names = FALSE,\n stringsAsFactors = FALSE\n)")
      
      thiscode <- c(
        "library(antitrust)",
        "\n\n ## Load Data:\n",
        indata_code,
        "\n\n ## Run Simulation: \n",
        atrfun,
        "\n\n ## Summary Tab Results:\n",
        "summary(simres, revenues = FALSE, levels = FALSE, market=TRUE)",
        "\n\n ## Details Tab Results:\n",
        "summary(simres, revenues = FALSE, levels = FALSE, market=FALSE)\n\n",
        "\n\n ## Elasticities Tab Results  (Pre-merger Only):\n",
        "elast(simres, preMerger = TRUE, market=TRUE)\n elast(simres, preMerger = TRUE, market=FALSE)",
        "\n\n ## Diagnostics Tab Results:\n",
        "calcDiagnostics(simres)\n\n"
      )
      
      cat(thiscode)
      
      
    })  
    
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

