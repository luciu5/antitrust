setClass(
  
  Class = "Cournot",
  contains="Bertrand",
  representation=representation(
    
    intercepts       = "vector",
    mcparm           = "vector",
    prices           = "vector",
    quantities       = "numeric",
    margins          = "numeric",
    quantityPre      = "numeric",
    quantityPost     = "numeric",
    productsPre      = "matrixOrVector",
    productsPost     = "matrixOrVector",
    demand           = "vector"
    
  ), 
  prototype(
    intercepts    =  numeric(),
    demand        =  "linear"
  ),
  validity=function(object){
    
    nprods <- length(object@shares) # count the number of products

    
    if(is.matrix(object@productsPre)){
      
      if(nprods != ncol(object@productsPre)){
        stop("The number of rows and columns in 'productsPre' must equal the length of 'labels'")}
      if(nrow(object@productsPre) != ncol(object@productsPre)){
        stop("'productsPre' must be a square matrix ")}
    }
    
    else if (nprods != length(object@productsPre)) stop("'productsPre' and 'labels' must be vectors of the same length")
    if(is.matrix(object@productsPost)){
      if(nprods != ncol(object@productsPost)){
        stop("The number of rows and columns in 'productsPost' must equal the length of 'labels'")}
      if(nrow(object@productsPost) != ncol(object@productsPost)){
        stop("'productsPost' must be a square matrix")}
    }
    
    else if (nprods != length(object@productsPost)) stop("'productsPost' and 'labels' must be vectors of the same length")
    
    
    
    if(nprods != length(object@quantities) ||
       nprods != length(object@margins) ||
       nprods != length(object@prices)){
      stop("'prices', 'quantities', 'margins', and 'shares' must all be vectors with the same length")}
    
    if(any(object@prices<0,na.rm=TRUE))             stop("'prices' values must be positive")
    if(any(object@quantities<0,na.rm=TRUE))          stop("'quantities' values must be positive")
    if(any(object@margins<0 | object@margins>1,na.rm=TRUE)) stop("'margins' values must be between 0 and 1")
    
    if(!tolower(trimws(object@demand)) %in% c("linear","log")){stop("'demand' must equal 'linear' or 'log'")}
    
    
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
    
    nprods <- length(products)
    quantTot <- tapply(quantities,products,sum)[products]
    quantOwn <- as.vector(owners %*% quantities)
    
    
    shares <- quantities/quantTot
    
    minDemand <- fun(theta){
      
      thisints <- theta[1:nprods]
      thisslopes <- theta[-(1:nprods)]
      
      elast <- ifelse(demand=="linear", prices*quantTot/thisslopes, thisslopes)
      
      FOC <- margins + shares/elast
      dem <- ifelse(demand=="linear", thisints + thisslopes*quantTot, 
                                     exp(thisints)*quantTot^thisslopes
                    )
      dem <- 1 - dem/prices
      dist <- c(FOC,dem)
      
      return(sum(dist^2,na.rm=TRUE))
    }
    
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
