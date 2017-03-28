setClass(
  
  Class = "Stackelberg",
  contains="Cournot",
  representation=representation(
    isLeaderPre = "matrix",
    isLeaderPost = "matrix",
    dmcfunPre           = "list",
    dmcfunPost          = "list"
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
      quantities <- object@quantityPre
      owner <- object@ownerPre}
    else{
      quantities <- object@quantityPost
      owner <- object@ownerPost}
    
    quantities[is.na(quantities)] <- 0
    
    quantOwner <- owner %*% quantities
    
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
      
      sharesOwner <- t( t(quantOwner) / mktQuant)
      
      
      elast <- 1/(t(quantOwner)/(prices/slopes)) * isLinear +
        (slopes / t(sharesOwner))  * ( 1 - isLinear)
      
      elast <- t(elast)
      
      dimnames(elast) <- object@labels
    }
    
    return(elast)
    
  }
)

## compute margins
setMethod(
  f= "calcMargins",
  signature= "Stackelberg",
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
  signature= "Stackelberg",
  definition=function(object){
    
    prices <- object@prices
    quantities <- object@quantities
    quantities[is.na(quantities)] <- 0
    margins <- object@margins
    products <- object@productsPre
    demand <- object@demand
    owner <- object@ownerPre
    vcfunPre <- object@vcfunPre
    nprods <- ncol(quantities)
    nplants <- nrow(quantities)
    
    noCosts <- length(vcfunPre) == 0
    isLeader <- object@isLeaderPre
    isLinear <- demand=="linear"
    
    quantTot <- colSums(quantities, na.rm = TRUE)
    quantPlants <- rowSums(quantities, na.rm = TRUE)
    
    quantOwner <- owner %*% quantities
    
    
    sharesOwner <- t(t(quantOwner)/quantTot)
    
    ## if variable cost function is missing
    ## assume vc_i = q_i^2/k_i for plant i
    ## adding create functions for mc and dmc
    if(!noCosts){
      mcPre <- sapply(1:nplants, function(i){object@mcfunPre(quantities[i,])})
      dmcPre <- sapply(1:nplants, function(i){object@dmcfunPre(quantities[i,])})
    }
    
    minDemand <- function(theta){
      
      if(noCosts){
        
        thiscap <- theta[1:nplants]
        theta <- theta[-(1:nplants)]
        mcPre <- quantPlants/thiscap
        dmcPre <- 1/thiscap
      }
      
      thisints <- theta[1:nprods]
      thisslopes <- theta[-(1:nprods)]
      
      thisprices <- ifelse(isLinear, thisints + thisslopes*quantTot, 
                           exp(thisints)*quantTot^thisslopes)
      
      thisPartial <- ifelse(isLinear, 
                            slopes,
                            exp(intercepts)*slopes*mktQuant^(slopes - 1))
      
      thisPass <- ifelse(isLeader,
                         sum(-thisPartial/(2*thisPartial - dmcPre)),
                         )
      thisPass[!isLeader] <- 0
      
      thisFOC <- (t(quantities) * thisPartial) %*% owner  -  mcPre
      thisFOC <- t(t(thisFOC) + thisprices + thisPartial*thisPass)
      
      
      
      
  
      
      
      
      dist <- c(FOC,prices)
      
      return(sum(dist^2,na.rm=TRUE))
    }
    
    
    bStart      =   ifelse(isLinear, -(prices*margins)/(sharesOwner*quantTot), -sharesOwner/margins)
    intStart    =   ifelse(isLinear,prices - bStart*quantTot, log(prices/(quantTot^bStart)))
    intStart    =   abs(intStart)
    
    parmStart   =   c( intStart,bStart)
    
    if(noCosts){
     capStart <-  
    }
    
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
    
    if(length(vcfunPre) ==0){
      
      elast <- 1/(t(quantOwner)/(prices/slopes)) * isLinear +
        (slopes / t(sharesOwner))  * ( 1 - isLinear)
      
      
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
    if(length(object@vcfunPost)==0){
      object@mcfunPost <- object@mcfunPre
      object@vcfunPost <- object@vcfunPre}
    
    object@intercepts <- intercepts
    object@slopes <-     slopes
    
    
    return(object)
    
  })    



setMethod(
  f= "calcQuantities",
  signature= "Stackelberg",
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




stackelberg <- function(prices,quantities,margins, 
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

  result <- new("Stackelberg",prices=prices, quantities=quantities,margins=margins,
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

