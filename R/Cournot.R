setClass(
  
  Class = "Cournot",
  contains="Bertrand",
  representation=representation(
    
    intercepts       = "vector",
    mcparm           = "vector",
    prices           = "vector",
    quantities       = "matrix",
    margins          = "matrix",
    quantityPre      = "matrix",
    quantityPost     = "matrix",
    quantityStart     = "numeric",
    productsPre      = "matrix",
    productsPost     = "matrix",
    demand           = "vector"
    
  ), 
  prototype(
    intercepts    =  numeric(),
    demand        =  "linear"
  ),
  validity=function(object){
    
    nprods <- ncol(object@quantities) # count the number of products
    nfirms  <- nrow(object@prices)     # count the number of firms
  
   if(!is.logical(object@productsPre)) stop("'productsPre' must be a logical matrix")
   if(!is.logical(object@productsPost)) stop("'productsPost' must be a logical matrix")
   
   if (!identical(dim(object@quantities), dim(object@margins))) stop("'margins' and 'quantities' must be matrices of the same dimension")
   if (!identical(dim(object@quantities), dim(object@productsPre))) stop("'productsPre' and 'quantities' must be matrices of the same dimension")
   if (!identical(dim(object@quantities), dim(object@productsPost))) stop("'productsPost' and 'quantities' must be matrices of the same dimension")
     
   
   if (nprods*nfirms != length(object@labels)) stop("'labels' length must equal the number of elements in 'quantities'")
    
    
    
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

## Generate a bunch of generic functions


setGeneric (
  name= "prodToMatrix",
  def=function(object,...){standardGeneric("prodToMatrix")}
)
setGeneric (
  name= "prodToVec",
  def=function(object,...){standardGeneric("prodToVec")}
)


## create ownership matrix
setMethod(
  f= "prodToMatrix",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){
    
    
    ## transform productsPre/productsPost vector into matrix, when applicable
    
    if(preMerger) {thisProd <- object@productsPre}
    else{         thisProd <- object@productsPost}
    
    
    
    if(is.vector(thisProd) || is.factor(thisProd)){
      
      nprod <- length(object@labels)
      products <- as.numeric(factor(thisProd))
      thisProd <- matrix(0,ncol=nprod,nrow=nprod)
      
      
      for( p in unique(products)){
        thisProd[products == p, products == p] = 1
      }
      
      
    }
    
    
    return(thisProd)
    
  }
)


## convert ownership matrix to vector
setMethod(
  f= "prodToVec",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){
    
    
    ## transform productsPre/productsPost matrix into an ownership vector
    if(preMerger) {thisProd <- object@productsPre}
    else{         thisProd <- object@productsPost}
    
    
    if(is.matrix(thisProd)){
      
      thisProd <- unique(thisProd)
      thisProd <- as.numeric(thisProd>=0.5) * (1: nrow(thisProd))
      thisProd <- apply(thisProd,2,max)
      
    }
    
    
    return(as.numeric(thisProd))
  }
  
)

setMethod(
  f= "calcShares",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE,revenue=FALSE){
    
    products <- object@products
    quantities <- calcQuantities(object,preMerger)
    
    if (revenue){
      if(preMerger){ prices <- object@pricePre}
      else{          prices <- object@pricePost}
      
      totrev <- tapply(prices*quantities, products, sum)[products]
      return(prices*quantities/totrev)
    }
    
    else{
      totquant <- tapply(quantities, products, sum)[products]
      return(quantities/totquant)}
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
      products <- object@productsPre
      owner <- objects@ownerPre
      quantities <- object@quantitiesPre}
    else{products <- object@productsPost
      owner <- objects@ownerPost
      quantities <- object@quantitiesPost}
    
    prices <- calcPrices(object,preMerger=preMerger)
    
    mktQuant <-  tapply(quantities,products)
    
    ##dQdP
    partial <- ifelse(demand=="linear", 
           slopes,
           exp(intercepts)*slopes*mktQuant^(slopes - 1))
    
    ##dPdQ
    partial <- 1/partial
   
    slopes    <- object@slopes
    
    if(market){
      
      elast <- partial*mktQuant/prices
    }
    
    else{
      
     prodmat <- model.matrix(~0+x,data.frame(x=products))
     
     elast <-  as.vector((quantities%*%(owner * prodmat))/partial)
     elast <- diag(1/elast)
     dimnames(elast,object@labels,object@labels)
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
  
  elast <- -1/diag(elast)
  names(elast) <- object@labels
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
    
    ownersVec <- ownerToVec(object,preMerger=TRUE)
    
    nprods <- length(prices)
    quantTot <- tapply(quantities,products,sum)
    quantOwn <- as.vector(owners %*% quantities)
    
    
    shares <- quantities/quantTot
    
    minDemand <- function(theta){
      
      thisints <- theta[1:nprods]
      thisslopes <- theta[-(1:nprods)]
      
      elast <- ifelse(demand=="linear", prices/(quantTot*thisslopes)[products], thisslopes)
      
      FOC <- margins + shares/elast
      dem <- ifelse(demand=="linear", thisints + thisslopes*quantTot, 
                                     exp(thisints)*quantTot^thisslopes
                    )
      dem <- 1 - dem/prices
      dist <- c(FOC,dem)
      
      return(sum(dist^2,na.rm=TRUE))
    }
    
    
    bStart      =   ifelse(demand=="linear", (prices*margins)/(shares*quantTot), -shares/margins)
    parmStart   =   c( bStart*quantTot + 1,bStart)
    
    
    
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
  
    elast <- ifelse(demand=="linear", prices/(quantTot*slopes)[products], slopes)
    elast <- shares/elast
    mc <- prices*(1+elast)
    
    mcparm <- tapply(quantOwn/mc,ownersVec,mean,na.rm=TRUE)[ownersVec]
    
    object@intercepts <- intercepts
    object@slopes <-     slopes
    object@mcparm <-     mcparm
      
  return(object)
    
  })    
setMethod(
  f= "calcMC",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){
    
    if(preMerger){ 
    owner  <- object@ownerPre
    quantity  <- object@quantityPre
   }
    else{          owner <-  object@ownerPost
    quantity  <- object@quantityPost
    }
    
    mcparm <- object@mcparm
    quantOwner <- as.vector(quantity %*% owner)
    mc <- quantOwner/mcparm
    names(mc) <- object@labels
    
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
                   mc        <- object@mcPre}
    else{          owner <-  object@ownerPost
                   products <-  object@productsPost
                   mc        <- object@mcPost}
    
    
    prodmat <- model.matrix(~0+x,data.frame(x=products))
    
    nprods <- length(products)
    
    FOC <- function(quantCand){
      
      if(preMerger){ object@quantityPre  <- quantCand}
      else{          object@quantityPost <- quantCand}
      
      thisPrice <- calcPrices(object, preMerger= preMerger)
      
      thisMC <- calcMC(object, preMerger= preMerger) 
      
      mktQuant <- tapply(quantCand,products)
      
      thisPartial <- ifelse(object@demand=="linear", 
                            slopes,
                         exp(intercepts)*object@slopes*mktQuant^(slopes - 1))
      
      thisPartial <- 1/thisPartial[products]
      
      thisFOC <- matrix(quantCand*thisPartial,ncol = nprods,nrow = nprods,byrow = TRUE)
      thisFOC <- thisFOC %*% (owner * prodmat) + thisPrice - thisMC
      
      return(sum(thisFOC^2))
    }
    
  })


setMethod(
  f= "calcPrices",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){
  
    if(preMerger){
      products <- object@productsPre
      quantities <- object@quantityPre
      }
    else{
      products <- object@productsPost
      quantities <- object@quantityPost
    }
    
    intercepts <- object@intercepts
    slopes     <- object@slopes
  

    mktQuant <- tapply(quantities, products, sum)
    mktQuant <- mktQuant[products]
    
    prices <- ifelse(object@demand == "linear",
                     intercepts + slopes * mktQuant,
                     exp(intercepts) * mktQuant^slopes
                     )
    
    
    return(prices)
      
  })  




cournot <- function(prices,quantities,margins, 
                    demand = rep("linear",length(prices)),
                    productsPre=!is.na(quantities), 
                    productsPost=productsPre, 
                    ownerPre,ownerPost,
                    mcDelta =rep(0,nrow(quantities)),
                    quantStart=as.vector(quantities),
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
    cname <- colnames(quantitites)
    }
    labels <- interaction(expand.grid(rname,cname),sep=":")
  }
  
  result <- new("Cournot",prices=prices, quantities=quantities,margins=margins,
                shares=shares,mcDelta=mcDelta, subset=subset,
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
  
  
  ## Calculate marginal cost
  result@mcPre <-  calcMC(result,TRUE)
  result@mcPost <- calcMC(result,FALSE)
  
  result@pricePre  <- calcPrices(result,TRUE,...)
  result@pricePost <- calcPrices(result,FALSE,subset=subset,...)
  
  
  return(result)
  
}

