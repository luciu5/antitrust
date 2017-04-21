setClass(
  Class   = "CournotCap",
  contains="Cournot",
  representation=representation(
      capacities           = "numeric"
    
  ),
  
  
  validity=function(object){
    
    
    
    
    
    nplants <- nrow(object@quantities)
    
    
    if(nplants != length(object@capacities)){
      stop("capacities' must be a vector whose length equals the number of rows in 'quantities'")}
    
    if(any(is.na(object@capacities) |
           !is.finite(object@capacities) |
           object@capacities<0 ,na.rm=TRUE)){stop("'capacities' values must be positive, finite numbers")}
    
    
    
    return(TRUE)
    
  })




setMethod(
  f= "calcSlopes",
  signature= "CournotCap",
  definition=function(object){
    
    ## Uncover Demand Coefficents
    
    prices <- object@prices
    quantities <- object@quantities
    quantities[is.na(quantities)] <- 0
    margins <- object@margins
    capacities <- object@capacities
    
    
    isConstrained <- capacities >= rowSums(quantities)
    isConstrained <- isConstrained * !is.na(margins) 
    
    mc <- t(t(1 - margins) * prices)
    
    products <- object@productsPre
    demand <- object@demand
    owner <- object@ownerPre
    mcfunPre <- object@mcfunPre
    nprods <- ncol(quantities)
    nplants <- nrow(quantities)
    
    noCosts <- length(mcfunPre) == 0
    isLinear <- demand=="linear"
    
    quantTot <- colSums(quantities, na.rm = TRUE)
    quantPlants <- rowSums(quantities, na.rm = TRUE)
    quantOwner <- owner %*% quantities
    
    
    if(!noCosts){
      mcPre <- sapply(1:nplants, function(i){object@mcfunPre[[i]](quantities[i,])})
    }
    
    sharesOwner <- t(t(quantOwner)/quantTot)
    
    minDemand <- function(theta){
      
      if(noCosts){
        
        mcPre <- theta[1:nplants]
        theta <- theta[-(1:nplants)] #constant marginal costs
        
      }
      
      thisints <- theta[1:nprods]
      thisslopes <- theta[-(1:nprods)]
      
      thisprices <- ifelse(isLinear, thisints + thisslopes*quantTot, 
                           exp(thisints)*quantTot^thisslopes)
      
      thisPartial <- ifelse(isLinear, 
                            thisslopes,
                            exp(thisints)*thisslopes*quantTot^(thisslopes - 1))
      
      
      thisFOC <- (t(quantities) * thisPartial ) %*% owner  + thisprices 
      thisFOC <- t(thisFOC) -   mcPre  
      
      ## eliminate a product FOC in plants with a binding capacity constraint 
      thisFOC[!isConstrained]
      
      
      dist <- c(thisFOC,thisprices - prices)
      
      if(noCosts){ dist <- c(dist, mcPre - mc)}
      
      return(sum(dist^2,na.rm=TRUE))
    }
    
    
    bStart      =   ifelse(isLinear,
                           colMeans(-(prices*margins)/(sharesOwner*quantTot),na.rm=TRUE), 
                           colMeans(-margins/sharesOwner,na.rm=TRUE))
    intStart    =   ifelse(isLinear,
                           prices - bStart*quantTot, 
                           log(prices/(quantTot^bStart)))
    intStart    =   abs(intStart)
    
    parmStart   =   c( intStart,bStart)
    
    if(noCosts){
      
      if(isLinear){margStart <- rowMeans(-(sharesOwner*quantTot)/(prices/bStart),na.rm=TRUE)}
      else{margStart <- rowMeans(-sharesOwner*bStart,na.rm=TRUE)} 
      
      mcStart  <- abs(prices*(margStart - 1)) 
      capStart <- quantPlants/mcStart  
      parmStart <- c(capStart,parmStart)
    }
    
    
    bestParms=optim(parmStart,minDemand)$par
    
    
    
    
    ## if no marginal cost functions are supplied 
    ## assume that plant i's marginal cost is
    ## q_i/k_i, where k_i is calculated from FOC
    
    if(noCosts){
      
      mcparm <- bestParms[1:nplants]
      bestParms <- bestParms[-(1:nplants)]
      
      
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
    
    intercepts = bestParms[1:nprods]
    slopes = bestParms[-(1:nprods)]
    
    
    object@intercepts <- intercepts
    object@slopes <-     slopes
    
    
    return(object)
    
  }
)
