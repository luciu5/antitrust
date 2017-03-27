setClass(
  
  Class = "Stackelberg",
  contains="Cournot",
  representation=representation(
    isLeaderPre = "matrix",
    isLeaderPost = "matrix",
    mcfunPre           = "list",
    mcfunPost           = "list"
  ), 
  prototype(
  ),
  validity=function(object){
    
    nplants <- nrow(object@quantities) # count the number of plants
    nprods  <- length(object@prices)     # count the number of products
  
    if(!is.logical(object@isLeaderPre)) stop("'leaderPre' must be a logical matrix")
    if(!is.logical(object@isLeaderPost)) stop("'leaderPost' must be a logical matrix")
    
    if(!identical(dim(object@quantities), dim(object@isLeaderPre))){stop("'isLeaderPre' must be a logical matrix whose dimensions must equal 'quantities' ")}
    if(!identical(dim(object@quantities), dim(object@isLeaderPost))){stop("'isLeaderPost' must be a logical matrix whose dimensions must equal 'quantities'")}
    
    
    
  }
)

setGeneric (
  name= "calcPass",
  def=function(object,...){standardGeneric("calcPass")}
)


##
## Stackelberg Methods
##


setMethod(
  f="calcPass",
  signature = "Stackelberg",
  definition=function(object, preMerger=TRUE){
  
    islinear <- object@demand=="linear"
    
    if(preMerger){isLeader <- object@isLeaderPre}
    else{isLeader <- object@isLeaderPost}
    
    pass <- ifelse( islinear, )
    
    
  }
  
)
setMethod(
  f= "elast",
  signature= "Stackelberg",
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
    owner <- object@ownerPre
    mcfunPre <- object@mcfunPre
    nprods <- ncol(quantities)
    nplants <- nrow(quantities)
    
    isLinear <- demand=="linear"
    
    quantTot <- colSums(quantities, na.rm = TRUE)
    quantPlants <- rowSums(quantities, na.rm = TRUE)
    
    
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
  
    
    
    ## if no marginal cost functions are supplied 
    ## assume that plant i's marginal cost is
    ## q_i/k_i, where k_i is calculated from FOC
    
    if(length(mcfunPre) ==0){
      
      elast <- 1/(t(quantities)/(prices/slopes)) * isLinear +
        (slopes / t(shares))  * ( 1 - isLinear)
      
      
      marg <- -1/t(elast)
      mcpred <- t(prices*(1-t(marg)))
      #mcact<- t(prices*(1-t(margins)))
      
      #mcact[is.na(mcact)] <- mcpred[is.na(mcact)]
      mcact = mcpred
      
      mcparm <- rowMeans(quantPlants/mcact,na.rm=TRUE)
      
      mcdef <- "function(q,mcparm = %f){ return(sum(q, na.rm=TRUE) * mcparm)}"
      mcdef <- sprintf(mcdef,1/mcparm)
      mcdef <- lapply(mcdef, function(x){eval(parse(text=x ))})
    
      object@mcfunPre <- mcdef
      names(object@mcfunPre) <- object@labels[[1]]
      
      vcdef <- "function(q,mcparm = %f){ return(sum(q, na.rm=TRUE)^2 * mcparm / 2)}"
      vcdef <- sprintf(vcdef,1/mcparm)
      vcdef <- lapply(vcdef, function(x){eval(parse(text=x ))})
      
      object@vcfunPre <- vcdef
      names(object@vcfunPre) <- object@labels[[1]]
      
    }
    if(length(object@mcfunPost)==0){
      object@mcfunPost <- object@mcfunPre
      object@vcfunPost <- object@vcfunPre}
    
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
    
   
   
    
    nplants <- nrow(quantity)
    
    mc <- rep(NA, nplants)
    
    for(f in 1:nplants){
      mc[f] <- mcfun[[f]](quantity[f,])
    }
    
    if(!preMerger){mc <- mc*(1 + object@mcDelta)}
    
    names(mc) <- object@labels[[1]]
    
    return(mc)
  })    


setMethod(
  f= "calcVC",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){
    
    if(preMerger){ 
      
      quantity  <- object@quantityPre
      vcfun <- object@vcfunPre
    }
    else{          
      
      quantity  <- object@quantityPost
      vcfun <- object@vcfunPost
    }
    
    
    
    
    nplants <- nrow(quantity)
    
    vc <- rep(NA, nplants)
    
    for(f in 1:nplants){
      vc[f] <- vcfun[[f]](quantity[f,])
    }
    
    if(!preMerger){vc <- vc*(1 + object@mcDelta)}
    
    names(vc) <- object@labels[[1]]
    
    return(vc)
  })   





## compute producer surplus
setMethod(
  f= "calcProducerSurplus",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){
    
    
    if( preMerger) {
      prices <- object@pricePre
      quantities <- object@quantityPre
      
    }
    else{prices <- object@pricePost
    quantities <- object@quantityPost
    
    }
    
    
    
    vc <- calcVC(object, preMerger= preMerger)
    
    ps <- colSums(prices*t(quantities), na.rm=TRUE) - vc
    names(ps) <- object@labels[[1]]
    
    return(ps)
  }
  
)


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
    
    nprods <- ncol(products)
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
      
      
      thisFOC <- (t(quantCand) * thisPartial) %*% owner - thisMC
      thisFOC <- t(t(thisFOC) + thisPrice)
      
      return(as.vector(thisFOC))
    }
    
    
    quantityStart <- sqrt(object@quantityStart[products]) #constrain positive
    
    ## Find price changes that set FOCs equal to 0
    minResult <- BBsolve( quantityStart,FOC, quiet=TRUE,control=object@control.equ,...)
    
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
  shares[is.na(shares)] <- 0 
  if(preMerger) {owner <- object@ownerPre}
  else{owner <- object@ownerPre}
  
  hhi <- colSums((owner %*% (shares*100))^2, na.rm =TRUE)
  
  return(hhi)
    
})


setMethod(
  f= "cmcr",
  signature= "Cournot",
  definition=function(object){
    
    owner <- object@ownerPre
    isParty <- rowSums( abs(object@ownerPost - object@ownerPre) ) > 0
    
    shares <- calcShares(object,preMerger=TRUE,revenue=FALSE)
    shares[is.na(shares)] <- 0 
    shares <- owner %*% shares
    shares <- shares[isParty,,drop=FALSE]
    mktElast <- elast(object, preMerger= TRUE,market=TRUE)

    cmcr <- -2 * apply(shares,2,prod) / (mktElast * colSums(shares) )
    
    return(cmcr * 100)
  })


setMethod(
  f= "summary",
  signature= "Cournot",
  definition=function(object,market=TRUE,revenue=FALSE,shares=FALSE,levels=FALSE,parameters=FALSE,digits=2,...){
    
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
                      plant=ids$plant,outPre=as.vector(t(outPre)), 
                      outPost = as.vector(t(outPost)))
    
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
      results <- subset(out, select = c(product,pricePre,pricePost,priceDelta,outPre,outPost,outDelta ))  
    }
    
    else{
      results <- subset(out, select = c(isParty,product,plant, pricePre,pricePost,priceDelta,outPre,outPost,outDelta ))
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


setMethod(
  f= "CV",
  signature= "Cournot",
  definition=function(object){
    
    demand <- object@demand
    slopes    <- object@slopes
    intercepts <- object@intercepts
    quantityPre <- colSums(object@quantityPre, na.rm=TRUE)
    quantityPost <- colSums(object@quantityPost, na.rm=TRUE)
    pricePre <- object@pricePre
    pricePost <- object@pricePost
    
    result <- ifelse(demand =="linear",
                 -.5*(pricePre - pricePost)*(quantityPre - quantityPost),
                 exp(intercepts)/(slopes + 1) * (quantityPre^(slopes+1) - quantityPost^(slopes+1)) -  (quantityPre - quantityPost)* pricePre
                 )
    
   
    names(result) <-  object@labels[[2]]
    return(result)
  })

setMethod(
  f= "calcPricesHypoMon",
  signature= "Cournot",
  definition=function(object,plantIndex){
    
    nhypoplants <- length(plantIndex)
    intercept <- object@intercepts
    nprods <- length(intercept)
    slopes <- object@slopes
    quantityPre <- as.vector(object@quantityPre)
    demand <- object@demand
    
    
    ## how to deal with multiple products?
   stop("A work in progress!!")
    
    calcMonopolySurplus <- function(quantCand){
      
      
      quantityPre[plantIndex] <- quantCand
      quantCand <- matrix(quantityCand,ncol=nprods)
      object@quantityPre <- quantCand
      mktQuant <- colSums(quantCand, na.rm = TRUE)
      
      priceCand <- ifelse(demand == "linear",
                      intercept + slopes * mktQuant,
                      exp(intercept)*mktQuant^slopes)
      
      vcCand <- calcVC(object, preMerger=TRUE)
      vcCand <- vcCand[plantIndex]
      
      revCand <-  colSums(priceCand*t(quantCand[plantIndex,]), na.rm=TRUE)
      
      
      surplus <- sum(revCand - vcCand, na.rm =TRUE)
      
      return(sum(surplus))
    }
    
    if( nhypoplants > 1){
     
      maxResult <- optim(object@quantityPre[prodIndex],
                         calcMonopolySurplus,
                         method="L-BFGS-B",
                         lower = rep(0,nhypoplants)
                               )
      
      quantitiesHM <- maxResult$par
    }
    
    
    else{
      
      upperB <- sum(quantityPre,na.rm=TRUE)
      maxResult <- optimize(calcMonopolySurplus,c(0, upperB),maximum = TRUE)
      quantitiesHM <- maxResult$maximum
    }
    
    quantityPre[plantIndex]
    names(pricesHM) <- object@labels[prodIndex]
    
    return(pricesHM)
    
    
  })



cournot <- function(prices,quantities,margins, 
                    demand = rep("linear",length(prices)),
                    mcfunPre=list(),
                    mcfunPost=mcfunPre,
                    vcfunPre=list(),
                    vcfunPost=vcfunPre,
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
                mcfunPre=mcfunPre, mcfunPost=mcfunPost,vcfunPre=vcfunPre, vcfunPost=vcfunPost,
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

