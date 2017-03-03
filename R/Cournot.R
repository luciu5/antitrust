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
    intercepts    =  numeric()
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
   if (nfirms != length(object@labels[[1]])) stop("'labels' length must be a list whose first element is a vector whose length equals the number of firms")
    
    if (nfirms != length(object@labels[[2]])) stop("'labels' length must be a list whose 2nd element is a vector whose length equals the number of products")
    
    
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
    
    if(preMerger) quantities <- object@quantitiesPre
    else{ quantities <- object@quantitiesPost}
    
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
    
    demand <- object@demand
    slopes <- object@slopes
    intercepts <- object@intercepts
    
  
    if(preMerger){
      owner <- objects@ownerPre
      quantities <- object@quantitiesPre}
    else{
      owner <- objects@ownerPost
      quantities <- object@quantitiesPost}
    
    prices <- calcPrices(object,preMerger=preMerger)
    
    mktQuant <-  colSums(quantities,na.rm = TRUE)
    
    ##dQdP
    partial <- ifelse(demand=="linear", 
           slopes,
           exp(intercepts)*slopes*mktQuant^(slopes - 1))
    
    ##dPdQ
    partial <- 1/partial
   
    
    if(market){
      
      elast <- partial*mktQuant/prices
    }
    
    else{
      
     shares <- t( t(quantities) / mktQuant)
     
     elast <- ifelse(demand=="linear", 
                     1/t( t(quantities) / (partial*prices) ),
                     t( slopes / t(shares)))
     
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
    owners <- object@ownePre
    mcfunPre <- object@mcfunPre
    nprods <- ncol(quantities)
    
    quantTot <- colSums(quantities, na.rm = TRUE)
    quantOwn <- rowSums(quantities, na.rm = TRUE)
    
    
    
    shares <- t(t(quantities)/quantTot)
    
    minDemand <- function(theta){
      
      thisints <- theta[1:nprods]
      thisslopes <- theta[-(1:nprods)]
      
      thisprices <- ifelse(demand=="linear", thisints + thisslopes*quantTot, 
                    exp(thisints)*quantTot^thisslopes)
      
      elast <- ifelse(demand=="linear", 
                      1/t(t(quantities)/(thisprices/thisslopes)),
                      t(thisslopes / t(shares)))
      
      
      FOC <- margins + 1/elast
      
                    
      dem <- 1 - thisprices/prices
      dist <- c(FOC,dem)
      
      return(sum(dist^2,na.rm=TRUE))
    }
    
    
    bStart      =   ifelse(demand=="linear", (prices*margins)/(shares*quantTot), -shares/margins)
    parmStart   =   c( -bStart*quantTot + 1,bStart)
    
    
    
    ## constrain diagonal elements so that D'b >=0
    ## constrain off-diagonal elements to be non-negative.
    
    #ui          =  diag(length(parmStart))
    #ui[1:nprods,1:nprods] = t(diversion)
    
    #ci = rep(0,length(parmStart)) 
    
    
    
    #bestParms=constrOptim(parmStart,minD,grad=NULL,ui=ui,ci=ci,
    #                      control=object@control.slopes)
    
    bestParms=optim(parmStart,minDemand)
    
    intercepts = bestParms$par[1:nprods]
    slopes =intercepts =  bestParms$par[-(1:nprods)]
  
    elast <- ifelse(demand=="linear", 
                    1/t(t(quantities)/(thisprices/slopes)),
                    t(slopes / t(shares)))
    
    marg <- -1/elast
    mc <- t(prices*(1-t(marg)))
    
    ## if no marginal cost functions are supplied 
    ## assume that firm i's marginal cost is
    ## q_i/k_i, where k_i is calculated from FOC
    
    if(length(mcfunPre) ==0){
      mcparm <- rowMeans(quantOwn/mc,na.rm=TRUE)
      
      fndef <- "function(q,mcparm = %f){ return(q/mcparm)}"
      fndef <- sprintf(fndef,mcparm)
      fndef <- lapply(fndef, function(x){eval(parse(text=x ))})
    
      object@mcfunPre <- fundef
      names(object@mcfunPre) <- object@labels[[1]]
      
    }
    if(length(mcfunPost)==0){object@mcfunPost <- object@mcfunPre}
    
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
  definition=function(object,preMerger=TRUE){
    
    slopes <- object@slopes
    intercept <- object@intercepts
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
      
      
      thisFOC <- t( t(quantCand) / thisPartial ) %*% owner - thisMC
      thisFOC <- t(thisFOC) + thisPrice
      
      return(sum(thisFOC^2))
    }
    
    
    quantityStart <- object@quantityStart[products]
    
    ## Find price changes that set FOCs equal to 0
    minResult <- BBsolve( quantityStart,FOC,quiet=TRUE,control=object@control.equ,...)
    
    if(minResult$convergence != 0){warning("'calcPrices' nonlinear solver may not have successfully converged. 'BBsolve' reports: '",minResult$message,"'")}
    
    
    if(isMax){
      
      hess <- genD(FOC,minResult$par) #compute the numerical approximation of the FOC hessian at optimium
      hess <- hess$D[,1:hess$p]
      hess <- hess * (owner[products,products]>0)   #0 terms not under the control of a common owner
      
      state <- ifelse(preMerger,"Pre-merger","Post-merger")
      
      if(any(eigen(hess)$values>0)){warning("Hessian of first-order conditions is not positive definite. ",state," price vector may not maximize profits. Consider rerunning 'calcPrices' using different starting values")}
    }
    
    
    quantEst        <- rep(NA, length(products))
    quantEst[products] <- minResult$par
    quantEst <- matrix(quantEst,ncol = nprods)
    
    dimnames(quantEst) <- object@labels
    
    
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
  
  shares <- quantities/sum(quantities)
  
  
  
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
                shares=shares,mcDelta=mcDelta, subset=subset, demand = demand,
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
  
  
  result@pricePre  <- calcPrices(result,TRUE,...)
  result@pricePost <- calcPrices(result,FALSE,subset=subset,...)
  
  
  return(result)
  
}

