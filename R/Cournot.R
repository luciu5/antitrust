setClass(
  
  Class = "Cournot",
  contains="Bertrand",
  representation=representation(
    
    intercepts       = "numeric",
    mcfunPre           = "list",
    mcfunPost           = "list",
    prices           = "vector",
    quantities       = "matrix",
    margins          = "matrix",
    quantityPre      = "matrix",
    quantityPost     = "matrix",
    quantityStart     = "numeric",
    productsPre      = "matrix",
    productsPost     = "matrix",
    demand           = "character"
    
  ), 
  prototype(
    intercepts    =  numeric(),
    mcfunPre =list(),
    mcfunPost=list()
  ),
  validity=function(object){
    
    nprods <- ncol(object@quantities) # count the number of products
    nfirms  <- nrow(object@prices)     # count the number of firms
  
   if(!is.list(object@mcfunPre) || 
      !length(object@mcfunPre) %in% c(nfirms,0)) {stop("'mcfunPre' must be a list of functions whose length equals the number of firms")}
   if(length(object@mcfunPre) >0 && any(sapply(object@mcfunPre,class) != "function"))
   {stop("'mcfunPre' must be a list of functions")}
    
    if(!is.list(object@mcfunPost) || 
       !length(object@mcfunPost) %in% c(nfirms,0)) {stop("'mcfunPost' must be a list of functions whose length equals the number of firms")}
    if(length(object@mcfunPost) >0 && any(sapply(object@mcfunPost,class) != "function"))
    {stop("'mcfunPost' must be a list of functions")}
    
   if(!is.logical(object@productsPre)) stop("'productsPre' must be a logical matrix")
   if(!is.logical(object@productsPost)) stop("'productsPost' must be a logical matrix")
   
   if (!identical(dim(object@quantities), dim(object@margins))) stop("'margins' and 'quantities' must be matrices of the same dimension")
   if (!identical(dim(object@quantities), dim(object@productsPre))) stop("'productsPre' and 'quantities' must be matrices of the same dimension")
   if (!identical(dim(object@quantities), dim(object@productsPost))) stop("'productsPost' and 'quantities' must be matrices of the same dimension")
     
  
   if(!is.list(object@labels)) stop("'labels' must be a list") 
   if (isTRUE(nfirms != length(object@labels[[1]]))) stop("'labels' length must be a list whose first element is a vector whose length equals the number of firms")
    
    if (isTRUE(nfirms != length(object@labels[[2]]))) stop("'labels' length must be a list whose 2nd element is a vector whose length equals the number of products")
    
    
    if(any(object@prices<=0,na.rm=TRUE))             stop("'prices' values must be positive")
    if(any(object@quantities<=0,na.rm=TRUE))          stop("'quantities' values must be positive")
    if(any(object@margins<0 | object@margins>1,na.rm=TRUE)) stop("'margins' values must be between 0 and 1")
    if(any(colSums(object@margins,na.rm=TRUE) == 0)) stop("at least one firm margin must be supplied for each product")
    
    if(all(!(object@demand %in% c("linear","log")))){stop("'demand' must equal 'linear' or 'log'")}
    if(length(object@demand) != nprods) stop("the length of 'demand' must equal the number of products")
    
    if(length(object@prices) != nprods) stop("the length of 'prices' must equal the number of products")
  }
)



##
## Cournot Methods
##


setMethod(
  f= "calcShares",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE,revenue=FALSE){
    
    if(preMerger) quantities <- object@quantityPre
    else{ quantities <- object@quantityPost}
    
    if (revenue){
      if(preMerger){ prices <- object@pricePre}
      else{          prices <- object@pricePost}
      
      totrev <- rowSums(prices*t(quantities), na.rm = TRUE)
      return(t(prices*t(quantities)/totrev))
    }
    
    else{
      totquant <- colSums(quantities)
      return(t(t(quantities)/totquant))}
  }
)

setMethod(
  f= "elast",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE,market=TRUE){
    
    isLinear <- object@demand =="linear"
    slopes <- object@slopes
    intercepts <- object@intercepts
    
  
    if(preMerger){
      quantities <- object@quantitiesPre}
    else{
      quantities <- object@quantitiesPost}
    
    prices <- calcPrices(object,preMerger=preMerger)
    
    mktQuant <-  colSums(quantities,na.rm = TRUE)
    
    ##dQdP
    partial <- ifelse(isLinear, 
           slopes,
           exp(intercepts)*slopes*mktQuant^(slopes - 1))
    
    ##dPdQ
    partial <- 1/partial
   
    
    if(market){
      
      elast <- partial*mktQuant/prices
    }
    
    else{
      
     shares <- t( t(quantities) / mktQuant)
     
     
     elast <- 1/(t(quantities)/(prices/slopes)) * isLinear +
       (slopes / t(shares))  * ( 1 - isLinear)
     
     elast <- t(elast)
     
     dimnames(elast) <- object@labels
    }
    
    return(elast)
    
  }
)

## compute margins
setMethod(
  f= "calcMargins",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){
    
    
  elast <- elast(object, preMerger = preMerger, market=FALSE)
  
  elast <- -1/elast
  dimnames(elast) <- object@labels
  return(elast)
  }
)


## Create a method to recover marginal cost using
## demand parameters and supplied prices


setMethod(
  f= "calcSlopes",
  signature= "Cournot",
  definition=function(object){
    
    prices <- object@prices
    quantities <- object@quantities
    margins <- object@margins
    products <- object@productsPre
    demand <- object@demand
    mcfunPre <- object@mcfunPre
    nprods <- ncol(quantities)
    nfirms <- nrow(quantities)
    
    isLinear <- demand=="linear"
    
    quantTot <- colSums(quantities, na.rm = TRUE)
    quantOwn <- rowSums(quantities, na.rm = TRUE)
    
    
    
    shares <- t(t(quantities)/quantTot)
    
    minDemand <- function(theta){
      
      thisints <- theta[1:nprods]
      thisslopes <- theta[-(1:nprods)]
      
      thisprices <- ifelse(isLinear, thisints + thisslopes*quantTot, 
                    exp(thisints)*quantTot^thisslopes)
      
      elast <- 1/(t(quantities)/(thisprices/thisslopes)) * isLinear +
       (thisslopes / t(shares))  * ( 1 - isLinear)
      
      
      
      FOC <- margins + 1/t(elast)
      
                    
      dem <- 1 - thisprices/prices
      dist <- c(FOC,dem)
      
      return(sum(dist^2,na.rm=TRUE))
    }
    
    
    bStart      =   ifelse(isLinear, -(prices*margins)/(shares*quantTot), -shares/margins)
    intStart    =   ifelse(isLinear,prices - bStart*quantTot, log(prices/(quantTot^bStart)))
    intStart    =   abs(intStart)
    
    parmStart   =   c( intStart,bStart)
    
    
    
    ## constrain diagonal elements so that D'b >=0
    ## constrain off-diagonal elements to be non-negative.
    
    #ui          =  -diag(length(parmStart))
    #for(i in 1:nprods){
    #  ui[i,i] <- 1
    #  ui[i,i+nprods] <- quantTot[i]
      
    #}
    #ui[1:nprods,1:nprods] = t(diversion)
    
    #ci = rep(0,length(parmStart)) 
    
    
    
    #bestParms=constrOptim(parmStart,minD,grad=NULL,ui=ui,ci=ci,
    #                      control=object@control.slopes)
    
    bestParms=optim(parmStart,minDemand)
    
    intercepts = bestParms$par[1:nprods]
    slopes = bestParms$par[-(1:nprods)]
  
    elast <- 1/(t(quantities)/(prices/slopes)) * isLinear +
      (slopes / t(shares))  * ( 1 - isLinear)
    
    
    marg <- -1/t(elast)
    mc <- t(prices*(1-t(marg)))
    
    ## if no marginal cost functions are supplied 
    ## assume that firm i's marginal cost is
    ## q_i/k_i, where k_i is calculated from FOC
    
    if(length(mcfunPre) ==0){
      mcparm <- rowMeans(quantOwn/mc,na.rm=TRUE)
      
      fndef <- "function(q,mcparm = %f){ return(q/mcparm)}"
      fndef <- sprintf(fndef,1/mcparm)
      fndef <- lapply(fndef, function(x){eval(parse(text=x ))})
    
      object@mcfunPre <- fndef
      names(object@mcfunPre) <- object@labels[[1]]
      
    }
    if(length(object@mcfunPost)==0){object@mcfunPost <- object@mcfunPre}
    
    object@intercepts <- intercepts
    object@slopes <-     slopes
    
      
  return(object)
    
  })    
setMethod(
  f= "calcMC",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){
    
    if(preMerger){ 
   
    quantity  <- object@quantityPre
    mcfun <- object@mcfunPre
   }
    else{          
    
    quantity  <- object@quantityPost
    mcfun <- object@mcfunPost
    }
    
   
    quantOwner <- rowSums(quantity, na.rm=TRUE)
    
    nfirms <- length(quantOwner)
    
    mc <- rep(NA, nfirms)
    
    for(f in 1:nfirms){
      mc[f] <- mcfun[[f]](quantOwner[f])
    }
    
    if(!preMerger){mc <- mc*(1 + object@mcDelta)}
    
    names(mc) <- object@labels[[1]]
    
    return(mc)
  })    

setMethod(
  f= "calcQuantities",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE,...){
    
    slopes <- object@slopes
    intercepts <- object@intercepts
    
    if(preMerger){ owner  <- object@ownerPre
                   products  <- object@productsPre
                   }
    else{          owner <-  object@ownerPost
                   products <-  object@productsPost
                   }
    
    nprods <- ncol(quantities)
    
    quantvec <- as.vector(quantities)
    products <- as.vector(products)
    
    FOC <- function(quantCand){
      
      quantCand <- quantCand^2 # constrain positive
      
      allquant <- rep(0,length(products))
      allquant[products] <- quantCand
      quantCand <- matrix(allquant,ncol=nprods)
      
      if(preMerger){ object@quantityPre  <- quantCand}
      else{          object@quantityPost <- quantCand}
      
      thisPrice <- calcPrices(object, preMerger= preMerger)
      
      thisMC <- calcMC(object, preMerger= preMerger) 
      
      mktQuant <- colSums(quantCand, na.rm=TRUE)
      
      thisPartial <- ifelse(object@demand=="linear", 
                            slopes,
                         exp(intercepts)*slopes*mktQuant^(slopes - 1))
      
      
      thisFOC <- (t(quantCand) / thisPartial) %*% owner - thisMC
      thisFOC <- t(thisFOC) + thisPrice
      
      return(as.vector(thisFOC))
    }
    
    
    quantityStart <- sqrt(object@quantityStart[products]) #constrain positive
    
    ## Find price changes that set FOCs equal to 0
    minResult <- BBsolve( quantityStart,FOC,quiet=TRUE,control=object@control.equ,...)
    
    if(minResult$convergence != 0){warning("'calcQuantities' nonlinear solver may not have successfully converged. 'BBsolve' reports: '",minResult$message,"'")}
    
    quantEst        <- rep(NA, length(products))
    quantEst[products] <- minResult$par^2
    quantEst <- matrix(quantEst,ncol = nprods)
    
    dimnames(quantEst) <- object@labels
    
    return(quantEst)
  })


setMethod(
  f= "calcPrices",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){
  
    if(preMerger){
          
      quantities <- object@quantityPre
      }
    else{
      quantities <- object@quantityPost
    }
    
    intercepts <- object@intercepts
    slopes     <- object@slopes
  

    mktQuant <- colSums(quantities, na.rm=TRUE)
    
    prices <- ifelse(object@demand == "linear",
                     intercepts + slopes * mktQuant,
                     exp(intercepts) * mktQuant^slopes
                     )
    
    
    names(prices) <- object@labels[[2]]
    return(prices)
      
  })  

## compute margins
setMethod(
  f= "calcMargins",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){
    
    if(preMerger) {
      prices <- object@pricePre
      products <- object@productsPre}
    else{prices <- object@pricePost
         products <- object@productsPost}
    
    mc     <- calcMC(object, preMerger = TRUE)
    
    margins <-  1 - t(t(products*mc)/prices)
    margins[!products] <- NA
    
    
    dimnames(margins) <- object@labels
    
    return(margins)
  }
  
)




setMethod(
  f= "summary",
  signature= "Cournot",
  definition=function(object,revenue=TRUE,shares=TRUE,levels=FALSE,parameters=FALSE,digits=2,...){
    
    curWidth <-  getOption("width")
    
    
    pricePre   <-  object@pricePre
    pricePost  <-  object@pricePost
    priceDelta <- calcPriceDelta(object,levels=levels)
    if(!levels) priceDelta <- priceDelta *100
    
    if(!shares){
      outPre  <-  object@quantityPre
      outPost <-  object@quantityPost
      
      if(revenue){
        outPre <- t(pricePre*t(outPre))
        outPost <- t(pricePost*t(outPost))
      }
      
      sumlabels=paste("quantity",c("Pre","Post"),sep="")
    }
    
    else{
      if(!shares){warning("'shares' equals FALSE but 'calcQuantities' not defined. Reporting shares instead of quantities")}
      
      outPre  <-  calcShares(object,preMerger=TRUE,revenue=revenue) * 100
      outPost <-  calcShares(object,preMerger=FALSE,revenue=revenue) * 100
      
      sumlabels=paste("shares",c("Pre","Post"),sep="")
    }
    
    mcDelta <- object@mcDelta * 100
    
    if(levels){outDelta <- outPost - outPre}
    else{outDelta <- (outPost/outPre - 1) * 100}
    
    
    isParty <- as.numeric(rowSums( abs(object@ownerPost - object@ownerPre))>0)
    isParty <- factor(isParty,levels=0:1,labels=c(" ","*"))
    
    results <- data.frame(pricePre=pricePre,pricePost=pricePost,
                          priceDelta=priceDelta,outputPre=outPre,
                          outputPost=outPost,outputDelta=outDelta)
    
    colnames(results)[colnames(results) %in% c("outputPre","outputPost")] <- sumlabels
    
    if(sum(abs(mcDelta))>0) results <- cbind(results,mcDelta=mcDelta)
    
    
    rownames(results) <- paste(isParty,object@labels)
    
    sharesPost <- calcShares(object,FALSE,revenue)
    
    cat("\nMerger simulation results under '",class(object),"' demand:\n\n",sep="")
    
    options("width"=100) # this width ensures that everything gets printed on the same line
    print(round(results,digits),digits=digits)
    options("width"=curWidth) #restore to current width
    
    cat("\n\tNotes: '*' indicates merging parties' products.\n ")
    if(levels){cat("\tDeltas are level changes.\n")}
    else{cat("\tDeltas are percent changes.\n")}
    if(revenue){cat("\tOutput is based on revenues.\n")}
    else{cat("\tOutput is based on units sold.\n")}
    
    results <- cbind(isParty, results)
    
    cat("\n\nShare-Weighted Price Change:",round(sum(sharesPost*priceDelta,na.rm=TRUE),digits),sep="\t")
    ##Only compute cmcr if cmcr method doesn't yield an error
    thisCMCR <- tryCatch(cmcr(object),error=function(e) FALSE)
    if(!is.logical(thisCMCR)){
      cat("\nShare-Weighted CMCR:",round(sum(cmcr(object)*sharesPost[isParty=="*"],na.rm=TRUE)/sum(sharesPost[isParty=="*"],na.rm=TRUE),digits),sep="\t")
    } 
    
    ##Only compute upp if prices are supplied
    thisUPP <- tryCatch(upp(object),error=function(e) FALSE)
    if(!is.logical(thisUPP)){
      cat("\nShare-Weighted Pricing Pressure:",round(sum(thisUPP*sharesPost[isParty=="*"],na.rm=TRUE)/sum(sharesPost[isParty=="*"],na.rm=TRUE),digits),sep="\t")}
    
    ##Only compute CV if prices  are supplied
    thisCV <- tryCatch(CV(object,...),error=function(e) FALSE)
    if(!is.logical(thisCV)){
      cat("\nCompensating Variation (CV):",round(thisCV,digits),sep="\t")}
    
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
    
    rownames(results) <- object@labels
    return(invisible(results))
    
  })



cournot <- function(prices,quantities,margins, 
                    demand = rep("linear",length(prices)),
                    mcfunPre=list(),
                    mcfunPost=list(),
                    productsPre=!is.na(quantities), 
                    productsPost=productsPre, 
                    ownerPre,ownerPost,
                    mcDelta =rep(0,nrow(quantities)),
                    quantityStart=as.vector(quantities),
                    control.slopes,
                    labels,
                   ...
){
  
  shares <- as.vector(quantities/sum(quantities))
  
  
  
  if(missing(labels)){
  if(is.null(dimnames(quantities))){ 
    rname <- paste0("O",1:nrow(quantities))
    cname <- paste0("P",1:ncol(quantities))
  }
    else{rname <- rownames(quantities)
    cname <- colnames(quantities)
    }
    labels <- list(rname,cname)
    
  }

  result <- new("Cournot",prices=prices, quantities=quantities,margins=margins,
                shares=shares,mcDelta=mcDelta, subset= rep(TRUE,length(shares)), demand = demand,
                mcfunPre=mcfunPre, mcfunPost=mcfunPost,
                ownerPre=ownerPre,productsPre=productsPre,productsPost=productsPost,
                ownerPost=ownerPost, quantityStart=quantityStart,labels=labels)
  
  
  if(!missing(control.slopes)){
    result@control.slopes <- control.slopes
  }
  
  
  ## Convert ownership vectors to ownership matrices
  result@ownerPre  <- ownerToMatrix(result,TRUE)
  result@ownerPost <- ownerToMatrix(result,FALSE)
  
  ## Calculate Demand Slope Coefficients and Intercepts
  result <- calcSlopes(result)
  
  
  result@quantityPre  <- calcQuantities(result, preMerger = TRUE,...)
  result@quantityPost <- calcQuantities(result,preMerger = FALSE,...)
  
  result@pricePre  <- calcPrices(result, preMerger = TRUE)
  result@pricePost <- calcPrices(result,preMerger = FALSE)
  
  return(result)
  
}

