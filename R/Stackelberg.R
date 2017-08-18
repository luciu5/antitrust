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
    isLeaderPre = matrix(),
    isLeaderPost = matrix(),
    dmcfunPre = list(),
    dmcfunPost=list()
  ),
  validity=function(object){
    
    nplants <- nrow(object@quantities) # count the number of plants
    nprods  <- length(object@prices)     # count the number of products
  
    if(!is.logical(object@isLeaderPre)) stop("'leaderPre' must be a logical matrix")
    if(!is.logical(object@isLeaderPost)) stop("'leaderPost' must be a logical matrix")
    
    if(!identical(dim(object@quantities), dim(object@isLeaderPre))){stop("'isLeaderPre' must be a logical matrix whose dimensions must equal 'quantities' ")}
    if(!identical(dim(object@quantities), dim(object@isLeaderPost))){stop("'isLeaderPost' must be a logical matrix whose dimensions must equal 'quantities'")}
    
    if( 
      !length(object@dmcfunPre) %in% c(nplants,0)) {stop("'dmcfunPre' must be a list of functions whose length equals the number of plants")}
    if(length(object@dmcfunPre) >0 && any(sapply(object@dmcfunPre,class) != "function"))
    {stop("'dmcfunPre' must be a list of functions")}
    
    if(
      !length(object@dmcfunPost) %in% c(nplants,0)) {stop("'dmcfunPost' must be a list of functions whose length equals the number of plants")}
    if(length(object@dmcfunPost) >0 && any(sapply(object@dmfunPost,class) != "function"))
    {stop("'dmcfunPost' must be a list of functions")}
    
    
  }
)


setGeneric (
  name= "calcdMC",
  def=function(object,...){standardGeneric("calcdMC")}
)


# setGeneric (
#   name= "calcPass",
#   def=function(object,...){standardGeneric("calcPass")}
# )


##
## Stackelberg Methods
##

# 
# setMethod(
#   f="calcPass",
#   signature = "Stackelberg",
#   definition=function(object, preMerger=TRUE){
#   
#     isLinearD <- object@demand=="linear"
#     
#     if(preMerger){isLeader <- object@isLeaderPre}
#     else{isLeader <- object@isLeaderPost}
#     
#     pass <- ifelse( isLinearD, )
#     
#     
#   }
#   
# )


setMethod(
  f= "calcdMC",
  signature= "Stackelberg",
  definition=function(object,preMerger=TRUE){
    
    if(preMerger){ 
      
      quantity  <- object@quantityPre
      dmcfun <- object@dmcfunPre
    }
    else{          
      
      quantity  <- object@quantityPost
      dmcfun <- object@dmcfunPost
    }
    
    
    
    
    nplants <- nrow(quantity)
    
    dmc <- rep(NA, nplants)
    
    for(f in 1:nplants){
      dmc[f] <- dmcfun[[f]](quantity[f,])
    }
    
    if(!preMerger){dmc <- dmc*(1 + object@mcDelta)}
    
    names(dmc) <- object@labels[[1]]
    
    return(dmc)
  })   




setMethod(
  f= "calcSlopes",
  signature= "Stackelberg",
  definition=function(object){
    
    
    prices <- object@prices
    quantities <- object@quantities
    quantities[is.na(quantities)] <- 0
    margins <- object@margins
    
    mc <- t(t(1 - margins) * prices)
    
    products <- object@productsPre
    demand <- object@demand
    owner <- object@ownerPre
    vcfunPre <- object@vcfunPre
    nprods <- ncol(quantities)
    nplants <- nrow(quantities)
    
    noCosts <- length(vcfunPre) == 0
    isLeader <- object@isLeaderPre
    isLinearD <- demand=="linear"
    isLinearC <- object@cost=="linear"
    
    quantTot <- colSums(quantities, na.rm = TRUE)
    quantPlants <- rowSums(quantities, na.rm = TRUE)
    quantOwner <- owner %*% quantities
    
    
    
    sharesOwner <- t(t(quantOwner)/quantTot)
    
    ## if variable cost function is missing
    ## assume vc_i = q_i^2/k_i for plant i
    ## adding create functions for mc and dmc
    if(!noCosts){
      mcPre <- sapply(1:nplants, function(i){object@mcfunPre[[i]](quantities[i,])})
      dmcPre <- sapply(1:nplants, function(i){object@dmcfunPre[[i]](quantities[i,])})
    }
    
    minParms <- function(theta){
      
      if(noCosts){
        
        thiscap <- theta[1:nplants]
        theta <- theta[-(1:nplants)]
        mcPre <- ifelse(isLinearC, quantPlants/thiscap, thiscap)
        dmcPre <- ifelse(isLinearC, 1/thiscap, 0)
      }
      
      thisints <- theta[1:nprods]
      thisslopes <- theta[-(1:nprods)]
      
      thisprices <- ifelse(isLinearD, thisints + thisslopes*quantTot, 
                           exp(thisints)*quantTot^thisslopes)
      
      thisPartial <- ifelse(isLinearD, 
                            thisslopes,
                            exp(thisints)*thisslopes*quantTot^(thisslopes - 1))
      
      dthisPartial <- ifelse(isLinearD,
                             0,
                             exp(thisints)*thisslopes*(thisslopes - 1)*quantTot^(thisslopes - 2))
      
    
      demPass <- dthisPartial * t(!isLeader * quantOwner)  
      thisPass <- -t((thisPartial + demPass)/
                   (2*thisPartial  + t(t(demPass) - dmcPre)))
                         
      thisPass[isLeader] <- 0
      
      thisFOC <- (t(quantities) * thisPartial + t(isLeader * quantities) * thisPartial * colSums(thisPass)) %*% owner  + thisprices 
      thisFOC <- t(thisFOC) -   mcPre  
      
  
      dist <- c(thisFOC,thisprices - prices)
      
      if(noCosts){ dist <- c(dist, mcPre - mc)}
      
      return(sum(dist^2,na.rm=TRUE))
    }
    
    
    bStart      =   ifelse(isLinearD,
                           colMeans(-(prices*margins)/(sharesOwner*quantTot),na.rm=TRUE), 
                           colMeans(-margins/sharesOwner,na.rm=TRUE))
    intStart    =   ifelse(isLinearD,
                           prices - bStart*quantTot, 
                           log(prices/(quantTot^bStart)))
    intStart    =   abs(intStart)
    
    parmStart   =   c( intStart,bStart)
    
    if(noCosts){
    
      if(isLinearD){margStart <- rowMeans(-(sharesOwner*quantTot)/(prices/bStart),na.rm=TRUE)}
      else{margStart <- rowMeans(-sharesOwner*bStart,na.rm=TRUE)} 
                       
     mcStart  <- abs(prices*(margStart - 1)) 
     capStart <- ifelse(isLinearC, quantPlants/mcStart, mcStart)   
     parmStart <- c(capStart,parmStart)
    }
    
    bestParms=optim(parmStart,minParms)$par
    
    ## if no variable cost functions are supplied 
    ## assume that plant i's varialbe cost is
    ## q_i^2/(2*k_i), where k_i is calculated from FOC
    
    if(noCosts){
      
      mcparm <- bestParms[1:nplants]
      bestParms <- bestParms[-(1:nplants)]
      
      
      dmcdef <- ifelse(isLinearC,"function(q,mcparm = %f){ val <-   1/mcparm; return(val)}",
                       "function(q,mcparm = %f){ val <-   0; return(val)}")
      dmcdef <- sprintf(dmcdef,mcparm)
      dmcdef <- lapply(dmcdef, function(x){eval(parse(text=x ))})
      
      object@dmcfunPre <- dmcdef
      names(object@dmcfunPre) <- object@labels[[1]]
      
      mcdef <- ifelse(isLinearC,"function(q,mcparm = %f){ val <- sum(q, na.rm=TRUE) / mcparm; return(val)}",
                      "function(q,mcparm = %f){ val <- mcparm; return(val)}")
      mcdef <- sprintf(mcdef,mcparm)
      mcdef <- lapply(mcdef, function(x){eval(parse(text=x ))})
      
      object@mcfunPre <- mcdef
      names(object@mcfunPre) <- object@labels[[1]]
      
      vcdef <- ifelse(isLinearC,"function(q,mcparm = %f){  val <-  sum(q, na.rm=TRUE)^2 / (mcparm * 2); return(val)}",
                      "function(q,mcparm = %f){  val <-  sum(q, na.rm=TRUE) * mcparm; return(val)}")
      vcdef <- sprintf(vcdef,mcparm)
      vcdef <- lapply(vcdef, function(x){eval(parse(text=x ))})
      
      object@vcfunPre <- vcdef
      names(object@vcfunPre) <- object@labels[[1]]
      
      
    }
      
    
    intercepts = bestParms[1:nprods]
    slopes = bestParms[-(1:nprods)]
    
    
    
   
    
    if(length(object@vcfunPost)==0){
      object@dmcfunPost <- object@dmcfunPre
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
    isLinearD <- object@demand=="linear"
    
    
    if(preMerger){ owner  <- object@ownerPre
                   products  <- object@productsPre
                   isLeader <- object@isLeaderPre
                   }
    else{          owner <-  object@ownerPost
                   products <-  object@productsPost
                   isLeader <- object@isLeaderPost
                   }
    
    nprods <- ncol(products)
    isProducts <- rowSums(products) > 0
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
      thisdMC <- calcdMC(object, preMerger= preMerger) 
      
      mktQuant <- colSums(quantCand, na.rm=TRUE)
      ownerQuant <- owner %*% quantCand
      
      thisPartial <- ifelse(isLinearD, 
                            slopes,
                         exp(intercepts)*slopes*mktQuant^(slopes - 1))
      
      dthisPartial <- ifelse(isLinearD,
                             0,
                             exp(intercepts)*slopes*(slopes - 1)*mktQuant^(slopes - 2))
      
      
      demPass <- dthisPartial * t(!isLeader * ownerQuant)  
      thisPass <- -t((thisPartial + demPass)/
        (2*thisPartial  + t(t(demPass) - thisdMC)))
      
    
      thisPass[isLeader | !isProducts] <- 0
      
      
      thisFOC <- (t(quantCand) * thisPartial  + t(isLeader * quantCand) * thisPartial*colSums(thisPass)) %*% owner + thisPrice 
      thisFOC <- t(thisFOC) - thisMC 
      
      thisFOC <- thisFOC[isProducts,]
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
                    cost   =   rep("linear",nrow(quantities)),
                    isLeaderPre = matrix(FALSE,ncol = ncol(quantities), nrow= nrow(quantities)),
                    isLeaderPost= isLeaderPre,
                    mcfunPre=list(),
                    mcfunPost=mcfunPre,
                    vcfunPre=list(),
                    vcfunPost=vcfunPre,
                    dmcfunPre=list(),
                    dmcfunPost=dmcfunPre,
                    capacitiesPre = rep(Inf,nrow(quantities)),
                    capacitiesPost = capacitiesPre,
                    productsPre=!is.na(quantities), 
                    productsPost=productsPre, 
                    ownerPre,ownerPost,
                    mcDelta =rep(0,nrow(quantities)),
                    quantityStart=as.vector(quantities),
                    control.slopes,
                    control.equ,
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
                shares=shares,mcDelta=mcDelta, subset= rep(TRUE,length(shares)), demand = demand, cost = cost,
                mcfunPre=mcfunPre, mcfunPost=mcfunPost,vcfunPre=vcfunPre, vcfunPost=vcfunPost,
                dmcfunPre=dmcfunPre, dmcfunPost=dmcfunPost, isLeaderPre = isLeaderPre, isLeaderPost = isLeaderPost,
                ownerPre=ownerPre,productsPre=productsPre,productsPost=productsPost,
                capacitiesPre=capacitiesPre,capacitiesPost=capacitiesPost,
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

