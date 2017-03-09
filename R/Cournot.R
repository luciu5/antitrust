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
    
    nfirms <- nrow(object@quantities) # count the number of firms
    nprods  <- length(object@prices)     # count the number of products
  
   if(!is.list(object@mcfunPre) || 
      !length(object@mcfunPre) %in% c(nfirms,0)) {stop("'mcfunPre' must be a list of functions whose length equals the number of firms")}
   if(length(object@mcfunPre) >0 && any(sapply(object@mcfunPre,class) != "function"))
   {stop("'mcfunPre' must be a list of functions")}
    
    if(length(object@mcfunPre) ==0 && any(rowSums(object@margins) == 0)){stop("When 'mcfunPre' is not supplied, at least one margin per firm must be supplied")}
    
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

    if(ncol(object@quantities) != nprods) stop("the number of columns in 'quantities' must equal the number of products")
    
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
      quantities <- object@quantityPre}
    else{
      quantities <- object@quantityPost}
    
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



## Method to compute HHI
setMethod(
  f= "hhi",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE, revenue = FALSE){
  
  shares <- calcShares(object,preMerger=preMerger,revenue=revenue)
  
  hhi <- colSums((shares*100)^2, na.rm =TRUE)
  
  return(hhi)
    
})


setMethod(
  f= "cmcr",
  signature= "Cournot",
  definition=function(object){
    
    isParty <- rowSums( abs(object@ownerPost - object@ownerPre) ) > 0
    
    shares <- calcShares(object,preMerger=TRUE,revenue=FALSE)
    shares <- shares[isParty,]
    mktElast <- elast(object, preMerger= TRUE,market=TRUE)

    cmcr <- -2 * apply(shares,2,prod) / (mktElast * colSums(shares) )
    
    return(cmcr * 100)
  })


setMethod(
  f= "summary",
  signature= "Cournot",
  definition=function(object,market=TRUE,revenue=TRUE,shares=FALSE,levels=FALSE,parameters=FALSE,digits=2,...){
    
    if(market){nfirms <- 1}
    else{ nfirms <- nrow(object@quantities) }
    
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
      ids <- data.frame(owner = 1 ,product= object@labels[[2]])  
    }
    
    
    else{
      
      ids <- expand.grid(owner=object@labels[[1]], product=object@labels[[2]])
    }
    
    
    out <- data.frame(product=ids$product, 
                      owner=ids$owner,outPre=as.vector(t(outPre)), 
                      outPost = as.vector(t(outPost)))
    
    if(market) {out$owner <- NULL}
    else{
      out$isParty <- as.numeric(rowSums( abs(object@ownerPost - object@ownerPre))>0)
      out$isParty <- factor(out$isParty,levels=0:1,labels=c(" ","*"))
    }
  
    mcDelta <- object@mcDelta * 100
    
    if(levels){out$outDelta <- out$outPost - out$outPre}
    else{out$outDelta <- (out$outPost/out$outPre - 1) * 100}
    
    out$pricePre <- rep(pricePre,each=nfirms)
    out$pricePost <- rep(pricePost,each=nfirms)
    out$priceDelta <- rep(priceDelta, each=nfirms)
    
    if(market){
      results <- subset(out, select = c(product,pricePre,pricePost,priceDelta,outPre,outPost,outDelta ))  
    }
    
    else{
      results <- subset(out, select = c(isParty,product,owner, pricePre,pricePost,priceDelta,outPre,outPost,outDelta ))
    }
    
    colnames(results)[colnames(results) %in% c("outputPre","outputPost")] <- sumlabels
    
    if(!market && sum(abs(mcDelta))>0) results <- cbind(results,mcDelta=mcDelta)
    
    
    
    sharesPost <- calcShares(object,FALSE,revenue)
    
    cat("\nMerger simulation results under '",class(object),"' demand:\n\n",sep="")
    
    options("width"=100) # this width ensures that everything gets printed on the same line
    options("scipen"=999) # this width ensures that everything gets printed on the same line
    print(cbind(results[,1:3], round(results[,-(1:3)], digits=digits)),digits=digits)
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
      
      print(round(cmcr(object),digits))
    } 
    
    ##Only compute upp if prices are supplied
    #thisUPP <- tryCatch(upp(object),error=function(e) FALSE)
    #if(!is.logical(thisUPP)){
    #  cat("\nShare-Weighted Pricing Pressure:",round(sum(thisUPP*sharesPost[isParty=="*"],na.rm=TRUE)/sum(sharesPost[isParty=="*"],na.rm=TRUE),digits),sep="\t")}
    
    ##Only compute CV if prices  are supplied
   # thisCV <- tryCatch(CV(object,...),error=function(e) FALSE)
   # if(!is.logical(thisCV)){
   #    cat("\nCompensating Variation (CV):",round(thisCV,digits),sep="\t")}
    
    cat("\n\n")
    
    
    if(parameters){
      
      cat("\nDemand Parameter Estimates:\n\n")
     
        print(round(object@slopes,digits))
        
      cat("\n\n")
      
      if(.hasSlot(object,"intercepts")){
        
        cat("\nIntercepts:\n\n")
        print(round(object@intercepts,digits))
        cat("\n\n")
        
      }
      
      
    }
    
    return(invisible(results))
    
  })


setMethod(
  f= "calcPricesHypoMon",
  signature= "Cournot",
  definition=function(object,firmIndex){
    
    nhypofirms <- length(firmIndex)
    intercept <- object@intercepts
    nprods <- length(intercept)
    slopes <- object@slopes
    mc <- object@mcPre[prodIndex]
    quantityPre <- as.vector(object@quantityPre)
    demand <- object@demand
    
    stop("This analysis is currently unavailable")
    
    
    
    calcMonopolySurplus <- function(quantCand){
      
      
      quantityPre[firmIndex] <- quantCand
      quantCand <- matrix(quantityCand,ncol=nprods)
      object@quantityPre <- quantCand
      mktQuant <- colSums(quantCand, na.rm = TRUE)
      
      priceCand <- ifelse(demand == "linear",
                      intercept + slopes * mktQuant,
                      exp(intercept)*mktQuant^slopes)
      
      mc <- calcMC(object, preMerger=TRUE)
      
      surplus <- sum(priceCand*t(quantityCand[firmIndex,]), na.rm =TRUE)
      
      return(sum(surplus))
    }
    
    ##Find starting value that always meets boundary conditions
    ##Note: if  nhypofirms=1, need to use a more accurate optimizer.
    
    if( nhypofirms> 1){
      
      if(det(slopes)!=0){startParm <- as.vector(solve(slopes) %*% (1 - intercept ))}
      else{startParm <- rep(0,nprods)}
      
      
      priceConstr <- pricePre
      priceConstr[prodIndex] <- 0
      
      maxResult <- constrOptim(startParm[prodIndex],calcMonopolySurplus,
                               grad=NULL,
                               ui=slopes[prodIndex,prodIndex],
                               ci=-intercept[prodIndex] - as.vector(slopes %*% priceConstr)[prodIndex],
                               control=list(fnscale=-1))
      
      pricesHM <- maxResult$par
    }
    
    
    else{
      
      upperB <- -(intercept[prodIndex] + sum(pricePre[-prodIndex]*slopes[prodIndex,-prodIndex]))/slopes[prodIndex,prodIndex]
      
      maxResult <- optimize(calcMonopolySurplus,c(0,upperB),maximum = TRUE)
      pricesHM <- maxResult$maximum
    }
    
    #priceDelta <- pricesHM/pricePre[prodIndex] - 1
    #names(priceDelta) <- object@labels[prodIndex]
    names(pricesHM) <- object@labels[prodIndex]
    
    return(pricesHM)
    
    
  })



cournot <- function(prices,quantities,margins, 
                    demand = rep("linear",length(prices)),
                    mcfunPre=list(),
                    mcfunPost=mcfunPre,
                    productsPre=!is.na(quantities), 
                    productsPost=productsPre, 
                    ownerPre,ownerPost,
                    mcDelta =rep(0,nrow(quantities)),
                    quantityStart=as.vector(quantities),
                    control.slopes,
                    control.eq,
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

