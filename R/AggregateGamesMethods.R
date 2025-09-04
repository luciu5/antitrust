#' @title Methods for Calculating Prices and Margins Using Aggregative Games
#' @name MarginsAG-Methods
#' @docType methods
#'
#' @aliases calcMarginsAG calcPricesAG
#' calcMarginsAG,ANY-method
#' calcPricesAG,ANY-method
#' calcMarginsAG,Logit-method
#' calcMarginsAG,CES-method
#' calcPricesAG,Logit-method
#'
#' @description Computes equilibrium product margins and prices using the aggregative games technique described in
#' Nocke and Schutz (2018). Assumes that firms are playing a
#' Nash-Bertrand pricing game with either Logit or CES demand
#' @param object An instance of one of the classes listed above.
#' @param preMerger If TRUE, returns pre-merger outcome. If
#' FALSE, returns post-merger outcome.  Default is TRUE.
#' @param level IF TRUE, return margins in dollars. If FALSE, returns 
#' margins in proportions. Default for most classes is FALSE.
#' @param subset A vector whose length equals the number  of products where each element equals TRUE if
#' the product indexed by that element should be included in the
#' post-merger simulation and FALSE if it should be excluded. Default is a
#' length k vector of TRUE.
#' @param isMax If TRUE, a check is run to determine if the calculated equilibrium price vector locally maximizes profits.
#' Default is FALSE.
#' @include MarginsMethods.R
#' @keywords methods
#' @references 
#' Nocke, V. and Schutz, N. (2018), Multiproduct-Firm Oligopoly: An Aggregative Games Approach. Econometrica, 86: 523-557.\doi{10.3982/ECTA14720}/
NULL

setGeneric (
  name= "calcMarginsAG",
  def=function(object,...){standardGeneric("calcMarginsAG")}
)

setGeneric (
  name= "calcPricesAG",
  def=function(object,...){standardGeneric("calcPricesAG")}
)

## compute margins using aggregate games method
## 
#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcMarginsAG",
  signature= "Logit",
  definition=function(object,preMerger=TRUE, level=FALSE){
    
    output <- ifelse(object@output,-1,1)
    meanval <- object@slopes$meanval
    alpha <- object@slopes$alpha
    idx <- object@normIndex
    
    outPrice <- object@priceOutside
    
    H0 <- is.na(idx)*exp(outPrice*alpha)
    
    if( preMerger) {
      mc <- object@mcPre
      owner  <- object@ownerPre
    }
    
    else{
      mc <- object@mcPost
      owner  <- object@ownerPost
    }
    
    
    type <-  exp(meanval + alpha*mc)
    
    typeFirm <- as.numeric(owner %*% type)
    
    mufun <- function(m,H){
      
      return(m*(1-(typeFirm/H)*exp(-m)) - 1)
    }
    
    Hfun <- function(h){
      
      muStart <- pmax(1,log(typeFirm/h))
      
      mu <- BB::BBsolve(muStart,mufun,H=h,quiet=TRUE)
    
      price <-  output*(mc - mu$par/alpha)
     
      val <- exp(meanval + alpha*price)
    
      omega <- sum(H0,val)/h
      
      return((omega - 1)^2)
    }
    
    HStart <- H0+sum(type)/exp(1)
    
    HBest <- optim(HStart,Hfun,method="Brent",lower=0,upper=1e6*HStart)
    
    
    muStart <- pmax(1,log(typeFirm/HBest$par))
    margins <- BB::BBsolve(muStart,mufun,H=HBest$par,quiet=TRUE)$par
    
    margins <- output*margins/alpha
    
    if(!level) {
      price <-  mc - margins
      margins <- margins / price }
    
    names(margins) <- object@labels
    
    return(as.numeric(margins))
  }
  
)

#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcMarginsAG",
  signature= "CES",
  definition=function(object,preMerger=TRUE, level=FALSE){
    
    output <- ifelse(object@output,1,-1)
    meanval <- object@slopes$meanval
    gamma <- object@slopes$gamma
    idx <- object@normIndex
    
    outPrice <- object@priceOutside
    
    H0 <- is.na(idx)*(outPrice^(1-gamma))
    
    if( preMerger) {
      mc <- object@mcPre
      owner  <- object@ownerPre
    }
    
    else{
      mc <- object@mcPost
      owner  <- object@ownerPost
    }
    
    
    type <-  meanval * mc^(1-gamma)
    
    typeFirm <- as.numeric(owner %*% type)
    
    mufun <- function(m,H){
      
      return(m*(1 - (gamma - 1)/gamma*(typeFirm/H)*(1 - m/gamma)^(gamma - 1)) - 1)
    }
    
    Hfun <- function(h){
      
      muStart <- pmax(1,gamma*(1 - (typeFirm/h)^(1/(gamma - 1))))
      
      mu <- BB::BBsolve(muStart,mufun,H=h,quiet=TRUE)
      
      price <-  mc/(1 - mu$par/gamma)
      
      val <- meanval * price^(1-gamma) 
      
      omega <- sum(H0,val)/h
      
      return((omega - 1)^2)
    }
    
    HStart <- H0+sum(type)*(1 - 1/gamma)^(gamma - 1)
    
    HBest <- optim(HStart,Hfun,method="Brent",lower=0,upper=1e6*HStart)
    
    
    muStart <- pmax(1,gamma*(1 - (typeFirm/HBest$par)^(1/(gamma - 1)))) 
    margins <- BB::BBsolve(muStart,mufun,H=HBest$par,quiet=TRUE)$par
    
    margins <- abs(margins/gamma)
    
    if(level) {
      price <- mc/(1 - output*margins)
      margins <- margins * price }
    
    names(margins) <- object@labels
    
    return(as.numeric(margins))
  }
  
)


## compute prices using aggregate games method
## 
#'@rdname MarginsAG-Methods
#'@export
setMethod(
  f= "calcPricesAG",
  signature= "Logit",
  definition=function(object,preMerger=TRUE,isMax=FALSE,subset){
    
    
    priceStart <- object@priceStart
    output    <-  object@output
    outSign <- ifelse(output,1,-1)
    
    if(preMerger){
      owner <- object@ownerPre
      mc    <- object@mcPre
    }
    else{
      owner <- object@ownerPost
      mc    <- object@mcPost
    }
    
  
      margins <- calcMarginsAG(object,preMerger=preMerger,level=TRUE)
      
      priceEst <- outSign*margins + mc
      
      names(priceEst) <- object@labels
      
      return(priceEst)
    }
)