#' @title Summary Methods
#' @description Summary methods for the \code{Bertrand}, \code{Auction2ndCap},  \code{Cournot}, and \code{Auction2ndLogit} classes.
#' Summarizes the effect of the merger, including price and revenue changes.
#' @name summary-methods
#' @docType methods
#'
#' @aliases summary,ANY-method
#' summary,AIDS-method
#' summary,Bertrand-method
#' summary,Auction2ndLogit-method
#' summary,Cournot-method
#' summary,Auction2ndCap-method
#' summary,VertBargBertLogit-method
#'
#' @param object an instance of class  \code{Bertrand}, \code{Auction2ndCap}, \code{Cournot}, or \code{Auction2ndLogit}
#' @param revenue When TRUE, returns revenues, when FALSE returns quantitities. Default is TRUE.
#' @param shares When TRUE, returns shares, when FALSE returns quantities (when possible). Default is TRUE.
#' @param levels When TRUE, returns changes in levels rather than percents and quantities rather than shares, when FALSE, returns
#' changes as a percent and shares rather than quantities. Default is FALSE.
#' @param parameters When TRUE, displays all demand parameters. Default is FALSE.
#' @param market When TRUE, displays aggregate information about the effect of a tariff.
#' When FALSE displays product-specific (or in the case of Cournot, plant-specific) effects.
#' Default is FALSE.
#' @param insideOnly When TRUE, rescales shares on inside goods to sum to 1. Default is FALSE.
#' @param digits Number of significant digits to report. Default is 2.
#' @param exAnte If \sQuote{exAnte} equals TRUE then the
#' \emph{ex ante} expected result for each firm is produced, while FALSE produces the
#' expected result conditional on each firm winning the auction. Default is FALSE.
#' @param ... Allows other objects to be passed to a \code{CV} method.
#'
#' @keywords methods
#' @include OwnershipMethods.R
NULL

#'@rdname summary-methods
#'@export
setMethod(
  f= "summary",
  signature= "Bertrand",
  definition=function(object,revenue=TRUE,shares=TRUE,levels=FALSE,parameters=FALSE,market=FALSE,insideOnly = TRUE,digits=2,...){

    curWidth <-  getOption("width")


    pricePre   <-  object@pricePre
    pricePost  <-  object@pricePost

    if(any(grepl("aids",class(object),ignore.case=TRUE))){

      priceDelta <-  object@priceDelta
    }
    else{ priceDelta <- calcPriceDelta(object,levels=levels)}

    if(!levels) priceDelta <- priceDelta *100

    if(!shares && !all(is.na(object@prices))){
      outPre  <-  calcQuantities(object,preMerger=TRUE)
      outPost <-  calcQuantities(object,preMerger=FALSE)

      if(revenue){
        outPre <- pricePre*outPre
        outPost <- pricePost*outPost
      }

      sumlabels=paste("quantity",c("Pre","Post"),sep="")
    }

    else{
      if(!shares){warning("'shares' equals FALSE but 'calcQuantities' not defined. Reporting shares instead of quantities")}

      outPre  <-  calcShares(object,preMerger=TRUE,revenue=revenue) * 100
      outPost <-  calcShares(object,preMerger=FALSE,revenue=revenue) * 100

      if(insideOnly){
        outPre <- outPre/sum(outPre)* 100
        outPost <- outPost/sum(outPost,na.rm=TRUE)* 100
      }

      sumlabels=paste("shares",c("Pre","Post"),sep="")
    }

    mcDelta <- object@mcDelta
    
    if(any(!grepl("auction2nd",class(object),ignore.case=TRUE))){
      
      mcDelta <- mcDelta * 100
    }
    

    if(levels){outDelta <- outPost - outPre}
    else{outDelta <- (outPost/outPre - 1) * 100}


    isParty <- as.numeric(rowSums( abs(object@ownerPost - object@ownerPre))>0)
    isParty <- factor(isParty,levels=0:1,labels=c(" ","*"))

    results <- data.frame(pricePre=pricePre,pricePost=pricePost,
                          priceDelta=priceDelta,outputPre=outPre,
                          outputPost=outPost,outputDelta=outDelta)



    if(sum(abs(mcDelta))>0) results <- cbind(results,mcDelta=mcDelta)


    rownames(results) <- paste(isParty,object@labels)

    sharesPost <- calcShares(object,FALSE,revenue)

    if(market){

      thiscmcr <- thiscv <- NA_real_

      try(thiscmcr <- cmcr(object,levels=levels), silent=TRUE)

      try(thiscv <- CV(object),silent = TRUE)

      thispsdelta  <- NA_real_
      try(thispsdelta  <- sum(calcProducerSurplus(object,preMerger=FALSE) - calcProducerSurplus(object,preMerger=TRUE),na.rm=TRUE),silent=TRUE)

      isparty <- isParty == "*"


      results <- with(results,
                      data.frame(
                        'HHI Change' = as.integer(HHI(outputPre/sum(outputPre),owner=object@ownerPost) - HHI(outputPre/sum(outputPre),owner=object@ownerPre)),
                        'Industry Price Change (%)' = sum(priceDelta * outputPost/sum(outputPost, na.rm = TRUE),na.rm=TRUE),
                        'Merging Party Price Change (%)'= sum(priceDelta[isparty] * outputPost[isparty], na.rm=TRUE) / sum(outputPost[isparty], na.rm=TRUE),
                        'Compensating Marginal Cost Reduction (%)' = sum(thiscmcr * outputPost[isparty]) / sum(outputPost[isparty], na.rm=TRUE),
                        'Consumer Harm ($)' = thiscv,
                        'Producer Benefit ($)' = thispsdelta,
                        'Difference ($)'= thiscv - thispsdelta,
                        check.names=FALSE
                      ))

      if(levels){colnames(results) <- gsub("%","$/unit",colnames(results))}


    }

    colnames(results)[colnames(results) %in% c("outputPre","outputPost")] <- sumlabels

    cat("\nMerger simulation results under '",class(object),"' demand:\n\n",sep="")

    options("width"=ifelse(market,25,100)) # this width ensures that everything gets printed on the same line
    print(round(results,digits),digits=digits, row.names=ifelse(market, FALSE, TRUE))
    options("width"=curWidth) #restore to current width



    if(!market){

      results <- cbind(isParty, results)
      rownames(results) <- object@labels

      cat("\n\tNotes: '*' indicates merging parties' products.\n ")
      if(levels){cat("\tDeltas are level changes.\n")}
      else{cat("\tDeltas are percent changes.\n")}
      if(revenue){cat("\tOutput is based on revenues.\n")}
      else{cat("\tOutput is based on units sold.\n")}

    }



    cat("\n\n")


    if(parameters){

      cat("\nDemand Parameter Estimates:\n\n")
      if(is.list(object@slopes)){
        print(lapply(object@slopes,round,digits=digits))
      }
      else{
        print(round(object@slopes,digits))
      }
      cat("\n\n")

      if(.hasSlot(object,"intercepts")){

        cat("\nIntercepts:\n\n")
        print(round(object@intercepts,digits))
        cat("\n\n")

      }

      if(.hasSlot(object,"constraint") && object@constraint){cat("\nNote: (non-singleton) nesting parameters are constrained to be equal")}
      cat("\n\n")

    }


    return(invisible(results))

  })

#'@rdname summary-methods
#'@export

setMethod(
  f= "summary",
  signature= "VertBargBertLogit",
  definition=function(object,revenue=TRUE,
                      levels=FALSE,parameters=FALSE,
                      market=FALSE,insideOnly = TRUE,
                      digits=2,...){
    
    curWidth <-  getOption("width")
    
    up <- object@up
    down <- object@down
    
    priceUpPre   <-  up@pricePre
    priceUpPost  <-  up@pricePost
    priceDownPre   <-  down@pricePre
    priceDownPost  <-  down@pricePost
    
    priceDelta <- calcPriceDelta(object,levels=levels,market=market)
    
    if(!levels) priceDelta <- lapply(priceDelta, function(x){x*100})
    
    #if(!shares && !all(is.na(object@prices))){
    #  outPre  <-  calcQuantities(object,preMerger=TRUE)
    #  outPost <-  calcQuantities(object,preMerger=FALSE)
      
    #  if(revenue){
    #    outPre <- pricePre*outPre
    #    outPost <- pricePost*outPost
    #  }
      
    #  sumlabels=paste("quantity",c("Pre","Post"),sep="")
    #}
    
    #else{
    #  if(!shares){warning("'shares' equals FALSE but 'calcQuantities' not defined.
    #Reporting shares instead of quantities")}
      
      outPre  <-  calcShares(object,preMerger=TRUE,revenue=revenue) * 100
      outPost <-  calcShares(object,preMerger=FALSE,revenue=revenue) * 100
      
      if(insideOnly){
        outPre <- outPre/sum(outPre)* 100
        outPost <- outPost/sum(outPost,na.rm=TRUE)* 100
      }
      
      sumlabels=paste("shares",c("Pre","Post"),sep="")
    #}
    
    mcDeltaUp   <- up@mcDelta * 100
    mcDeltaDown <- down@mcDelta * 100
    
    if(levels){outDelta <- outPost - outPre}
    else{outDelta <- (outPost/outPre - 1) * 100}
    
    
    isPartyHorzDown <- down@ownerPost %in% down@ownerPre &
                       down@ownerPost != down@ownerPre
    if(any(isPartyHorzDown)){
    isPartyHorzDown <- down@ownerPost==down@ownerPost[isPartyHorzDown]
}
isPartyHorzUp <- up@ownerPost %in% up@ownerPre &
                 up@ownerPost !=   up@ownerPre
if(any(isPartyHorzUp)){
  isPartyHorzUp <- up@ownerPost==up@ownerPost[isPartyHorzUp]
}
    isPartyVert <- unique(down@ownerPost[down@ownerPost == up@ownerPost &
                       down@ownerPre != up@ownerPre])
    isPartyVert <- (down@ownerPost %in% isPartyVert) | (up@ownerPost %in% isPartyVert)
  
    isParty <- factor(isPartyHorzDown | isPartyHorzUp |isPartyVert,levels=c(FALSE,TRUE),labels=c(" ","*"))
    
    results <- data.frame(priceUpPre=priceUpPre,
                          priceUpPost=priceUpPost,
                          priceUpDelta=priceDelta$up,
                          priceDownPre=priceDownPre,
                          priceDownPost=priceDownPost,
                          priceDownDelta=priceDelta$down,
                          outputPre=outPre,
                          outputPost=outPost,outputDelta=outDelta)
    
    
    
    if(sum(abs(mcDeltaUp))>0) results <- cbind(results,mcDeltaUp=mcDeltaUp)
    if(sum(abs(mcDeltaDown))>0) results <- cbind(results,mcDeltaDown=mcDeltaDown)
    
    
    rownames(results) <- paste(isParty,down@labels)
    
    sharesPost <- calcShares(object,FALSE,revenue)
    
    if(market){
      
      thiscmcr <- thiscv <- NA_real_
      #try(thiscmcr <- cmcr(object), silent=TRUE)
      try(thiscv <- CV(object),silent = TRUE)
      
      try(thispsPre <- calcProducerSurplus(object,TRUE),silent=TRUE)
      try(thispsPost <- calcProducerSurplus(object,FALSE),silent=TRUE)
      thispsdeltaUp  <- sum(thispsPost$up - thispsPre$up,na.rm=TRUE)
      thispsdeltaDown  <- sum(thispsPost$down - thispsPre$down,na.rm=TRUE)
      
      isparty <- isParty == "*"
      
      
      hhiUp <- as.integer(HHI(outPre/sum(outPre),owner=up@ownerPost) - HHI(outPre/sum(outPre),owner=up@ownerPre))
      hhiDown <- as.integer(HHI(outPre/sum(outPre),owner=down@ownerPost) - HHI(outPre/sum(outPre),owner=down@ownerPre))
      
      partylabel <- unique(down@ownerPost[isParty])
      downOwnerPost <- down@ownerPost
      downOwnerPost[isparty]=partylabel[1]
      hhiVert <- as.integer(HHI(outPre/sum(outPre),owner=downOwnerPost) - HHI(outPre/sum(outPre),owner=down@ownerPre))
      
      hhidelta <- ifelse(any(isPartyVert),hhiVert,max(hhiUp,hhiDown))
      
      results <- with(results,
                      data.frame(
                        'HHI Change' =  hhidelta,
                        'Up Price Change (%)' = priceDelta$up,
                        'Down Price Change (%)' = priceDelta$down,
                        #'Merging Party Price Change (%)'= sum(priceDelta[isparty] * outputPost[isparty], na.rm=TRUE) / sum(outputPost[isparty]),
                        #'Compensating Marginal Cost Reduction (%)' = sum(thiscmcr * outputPost[isparty]) / sum(outputPost[isparty]),
                        'Consumer Harm ($)' = thiscv,
                        'Up Producer Benefit ($)' = thispsdeltaUp,
                        'Down Producer Benefit ($)' = thispsdeltaDown,
                        'Difference ($)'= thiscv - thispsdeltaUp - thispsdeltaDown,
                        check.names=FALSE
                      ))
      
      if(levels){colnames(results) <- gsub("%","$/unit",colnames(results))}
      
      
    }
    
    colnames(results)[colnames(results) %in% c("outputPre","outputPost")] <- sumlabels
    
    cat("\nMerger simulation results under '",class(object),"' demand:\n\n",sep="")
    
    options("width"=ifelse(market,25,100)) # this width ensures that everything gets printed on the same line
    print(round(results,digits),digits=digits, row.names=ifelse(market, FALSE, TRUE))
    options("width"=curWidth) #restore to current width
    
    
    
    if(!market){
      
      results <- cbind(isParty, results)
      rownames(results) <- down@labels
      
      cat("\n\tNotes: '*' indicates merging parties' products.\n ")
      if(levels){cat("\tDeltas are level changes.\n")}
      else{cat("\tDeltas are percent changes.\n")}
      if(revenue){cat("\tOutput is based on revenues.\n")}
      else{cat("\tOutput is based on units sold.\n")}
      
    }
    
    
    
    cat("\n\n")
    
    
    if(parameters){
      
      print(getParms(object), digits=digits)
        
      
      if(.hasSlot(object,"constraint") && object@constraint){cat("\nNote: (non-singleton) nesting parameters are constrained to be equal")}
      cat("\n\n")
      
    }
    
    
    return(invisible(results))
    
  })



#'@rdname summary-methods
#'@export
setMethod(
  f= "summary",
  signature= "Auction2ndCap",
  definition=function(object,exAnte=FALSE,parameters=FALSE,market=TRUE,digits=2){

    curWidth <-  getOption("width")


    pricePre   <-  calcPrices(object,preMerger=TRUE,exAnte=exAnte)
    pricePost  <-  calcPrices(object,preMerger=FALSE,exAnte=exAnte)
    priceDelta <- (pricePost/pricePre - 1) * 100


    outPre  <-  calcShares(object,TRUE,exAnte=exAnte) * 100
    outPost <-  calcShares(object,FALSE,exAnte=exAnte) * 100





    mcDelta <- object@mcDelta

    outDelta <- (outPost/outPre - 1) * 100


    isParty <- object@ownerPost != object@ownerPre
    isParty <- c(object@ownerPre[isParty],object@ownerPost[isParty])
    isParty <- factor(ifelse(object@ownerPre %in% isParty,1,0),levels=0:1,labels=c(" ","*"))

    results <- data.frame(pricePre=pricePre,pricePost=pricePost,
                          priceDelta=priceDelta,sharesPre=outPre,
                          sharesPost=outPost,sharesDelta=outDelta)


    if(sum(abs(mcDelta))>0) results <- cbind(results,mcDelta=mcDelta)


    rownames(results) <- paste(isParty,object@labels)


    if( market){
    
      thiscmcr <- thiscv <- NA_real_
      #try(thiscmcr <- cmcr(object), silent=TRUE)
      try(thiscv <- CV(object),silent = TRUE)
      
      try(thispsPre <- calcProducerSurplus(object,TRUE),silent=TRUE)
      try(thispsPost <- calcProducerSurplus(object,FALSE),silent=TRUE)
      thispsdelta  <- sum(thispsPost - thispsPre,na.rm=TRUE)
      
      
      results <- with(results,
                    data.frame(
                      'HHI Change' = as.integer(HHI(outPre/sum(outPre),owner=object@ownerPost) - HHI(outPre/sum(outPre),owner=object@ownerPre)),
                      'Industry Price Change (%)' = sum(priceDelta * outPost/sum(outPost, na.rm = TRUE),na.rm=TRUE),
                      'Merging Party Price Change (%)'= sum(priceDelta[isParty] * outPost[isParty], na.rm=TRUE) / sum(outPost[isParty], na.rm=TRUE),
                      'Compensating Marginal Cost Reduction (%)' = sum(thiscmcr * outPost[isParty]) / sum(outPost[isParty], na.rm=TRUE),
                      'Consumer Harm ($)' = thiscv,
                      'Producer Benefit ($)' = thispsdelta,
                      'Difference ($)'= thiscv - thispsdelta,
                      check.names=FALSE
                    ))
    }
    cat("\nMerger simulation results under '",class(object),"':\n\n",sep="")

    options("width"=100) # this width ensures that everything gets printed on the same line
    print(round(results,digits),digits=digits)
    options("width"=curWidth) #restore to current width

    if(!market){
    cat("\n\tNotes: '*' indicates merging parties. Deltas are percent changes.\n")

    if(exAnte){cat("\tEx Ante shares and prices are reported.\n")}
    else{cat("\tShares and prices conditional on a firm winning are reported.\n")}

    results <- cbind(isParty, results)

    cat("\n\nPre-Merger Buyer Reserve:",round(object@reservePre,digits),sep="\t")
    cat("\nPost-Merger Buyer Reserve:",round(object@reservePost,digits),sep="\t")
    cat("\n\n% Change In Expected Price:",round((calcExpectedPrice(object,FALSE)-calcExpectedPrice(object,TRUE))/calcExpectedPrice(object,TRUE)*100,digits),sep="\t")
    cat("\n")
    cat("% Change In Buyer's Expected Cost:",round((calcBuyerExpectedCost(object,FALSE)-calcBuyerExpectedCost(object,TRUE))/calcBuyerExpectedCost(object,TRUE)*100,digits),sep="\t")
    cat("\n\n")

    rownames(results) <- object@labels
    }
    
    if(parameters){

      cat("\nSupplier Cost Distribution Parameters:\n\n")

      print(round(object@sellerCostParms,digits))

      cat("\nBuyer Valuation:\n\n")
      print(round(object@buyerValuation,digits))

      cat("\n\n")

    }

   
    return(invisible(results))

  })

#'@rdname summary-methods
#'@export
setMethod(
f= "summary",
signature= "Cournot",
definition=function(object,market=FALSE,revenue=FALSE,shares=FALSE,levels=FALSE,parameters=FALSE,digits=2,...){

  if(market){nplants <- 1}
  else{ nplants <- nrow(object@quantities) }

  curWidth <-  getOption("width")
  curSci  <-  getOption("scipen")

  pricePre   <-  object@pricePre
  pricePost  <-  object@pricePost
  priceDelta <- calcPriceDelta(object,levels=levels)
  if(!levels) priceDelta <- priceDelta *100

  if(!shares){
    outPre  <-  object@quantityPre
    outPost <-  object@quantityPost
    sumlabels=paste("quantity",c("Pre","Post"),sep="")

    if(revenue){
      outPre <- t(pricePre*t(outPre))
      outPost <- t(pricePost*t(outPost))
      sumlabels=paste("revenue",c("Pre","Post"),sep="")
    }

  }

  else{
    if(!shares){warning("'shares' equals FALSE but 'calcQuantities' not defined. Reporting shares instead of quantities")}


    outPre  <-  calcShares(object,preMerger=TRUE,revenue=revenue) * 100
    outPost <-  calcShares(object,preMerger=FALSE,revenue=revenue) * 100


    sumlabels=paste("shares",c("Pre","Post"),sep="")
  }

  if(market){

    outPre <- colSums(outPre,na.rm=TRUE)
    outPost <- colSums(outPost,na.rm=TRUE)
    ids <- data.frame(plant = 1 ,product= object@labels[[2]])
  }


  else{

    ids <- expand.grid(plant=object@labels[[1]], product=object@labels[[2]])
  }


  out <- data.frame(product=ids$product,
                    plant=ids$plant,outPre=as.vector(outPre),
                    outPost = as.vector(outPost))

  if(market) {out$plant <- NULL}
  else{
    out$isParty <- as.numeric(rowSums( abs(object@ownerPost - object@ownerPre))>0)
    out$isParty <- factor(out$isParty,levels=0:1,labels=c(" ","*"))
  }

  mcDelta <- object@mcDelta * 100

  if(levels){out$outDelta <- out$outPost - out$outPre}
  else{out$outDelta <- (out$outPost/out$outPre - 1) * 100}

  out$pricePre <- rep(pricePre,each=nplants)
  out$pricePost <- rep(pricePost,each=nplants)
  out$priceDelta <- rep(priceDelta, each=nplants)

  if(market){
    results <- out[,c("product","pricePre","pricePost","priceDelta","outPre","outPost","outDelta" )]
  }

  else{
    results <- out[, c("isParty","product","plant", "pricePre","pricePost","priceDelta","outPre","outPost","outDelta" )]
  }

  colnames(results)[colnames(results) %in% c("outPre","outPost")] <- sumlabels

  if(!market && sum(abs(mcDelta))>0) results <- cbind(results,mcDelta=mcDelta)



  sharesPost <- calcShares(object,FALSE,revenue)

  cat("\nMerger simulation results under '",class(object),"' demand:\n\n",sep="")

  options("width"=100) # this width ensures that everything gets printed on the same line
  options("scipen"=999) # this width ensures that everything gets printed on the same line
  print(format(results,digits=digits),row.names = FALSE)
  options("width"=curWidth) #restore to current width
  options("scipen"=curSci) #restore to current scientific notation

  cat("\n\tNotes: '*' indicates merging parties' products.\n ")
  if(levels){cat("\tDeltas are level changes.\n")}
  else{cat("\tDeltas are percent changes.\n")}
  if(revenue){cat("\tOutput is based on revenues.\n")}
  else{cat("\tOutput is based on units sold.\n")}



  ##Only compute cmcr if cmcr method doesn't yield an error
  thisCMCR <- tryCatch(cmcr(object),error=function(e) FALSE)
  if(!is.logical(thisCMCR)){
    cat("\n\nCMCR:\n\n")

    cat(format(cmcr(object),digits=digits), fill=TRUE,labels=object@labels[[2]])
  }

  ##Only compute upp if prices are supplied
  #thisUPP <- tryCatch(upp(object),error=function(e) FALSE)
  #if(!is.logical(thisUPP)){
  #  cat("\nShare-Weighted Pricing Pressure:",format(sum(thisUPP*sharesPost[isParty=="*"],na.rm=TRUE)/sum(sharesPost[isParty=="*"],na.rm=TRUE),digits),sep="\t")}

  ##Only compute CV if prices  are supplied
  thisCV <- tryCatch(CV(object,...),error=function(e) FALSE)
  if(!is.logical(thisCV)){
    cat("\n\nCompensating Variation (CV):\n\n")
    cat(format(thisCV,digits=digits),fill=TRUE, labels=object@labels[[2]])}

  cat("\n\n")


  if(parameters){

    cat("\nDemand Parameter Estimates:\n\n")

    print(format(object@slopes,digits=digits))

    cat("\n\n")

    if(.hasSlot(object,"intercepts")){

      cat("\nIntercepts:\n\n")
      print(format(object@intercepts,digits=digits))
      cat("\n\n")

    }


  }

  return(invisible(results))

})

#'@rdname summary-methods
#'@export
setMethod(
  f= "summary",
  signature= "Auction2ndLogit",
  definition=function(object,levels=TRUE,revenue=FALSE,...){

    callNextMethod(object,levels=levels,revenue=revenue,...)
  }

)
