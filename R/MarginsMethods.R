#' @title Methods for Calculating Diagnostics
#' @name Margins-Methods
#' @docType methods
#'
#' @aliases calcMargins
#' calcMargins,ANY-method
#' calcMargins,AIDS-method
#' calcMargins,Bertrand-method
#' calcMargins,Bargaining-method
#' calcMargins,LogitCap-method
#' calcMargins,Auction2ndLogit-method
#' calcMargins,Cournot-method
#'
#' @description Computes equilibrium product margins assuming that firms are playing a
#' Nash-Bertrand or Cournot game. For "LogitCap", assumes firms are
#' playing a Nash-Bertrand or Cournot game with capacity constraints.
#'
#' @param object An instance of one of the classes listed above.
#' @param preMerger If TRUE, returns pre-merger outcome. If
#' FALSE, returns post-merger outcome.  Default is TRUE.
#' @param exAnte If \sQuote{exAnte} equals TRUE then the
#' \emph{ex ante} expected result for each firm is produced, while FALSE produces the
#' expected result conditional on each firm winning the auction. Default is FALSE.
#' @param subset A vector of length k where each element equals TRUE if
#' the product indexed by that element should be included in the
#' post-merger simulation and FALSE if it should be excluded.Default is a
#' length k vector of TRUE.
#'
#' @include CostMethods.R
#' @keywords methods
NULL

setGeneric (
  name= "calcMargins",
  def=function(object,...){standardGeneric("calcMargins")}
)

## compute margins
#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "Bertrand",
  definition=function(object,preMerger=TRUE){



    if( preMerger) {

      owner  <- object@ownerPre
      revenue<- calcShares(object,preMerger,revenue=TRUE)

      elast <-  elast(object,preMerger)
      margins <-  -1 * as.vector(MASS::ginv(t(elast)*owner) %*% (revenue * diag(owner))) / revenue


    }

    else{
      prices <- object@pricePost
      mc     <- object@mcPost

      margins <- 1 - mc/prices
    }


    names(margins) <- object@labels

    return(as.vector(margins))
  }

)

#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "Bargaining",
  definition=function(object,preMerger=TRUE){
    
    
    
    if( preMerger) {
      
      owner  <- object@ownerPre
      revenue<- calcShares(object,preMerger,revenue=TRUE)
      
      elast <-  elast(object,preMerger)
      margins <-  -1 * as.vector(MASS::ginv(t(elast)*owner) %*% (revenue * diag(owner))) / revenue
      
      
    }
    
    else{
      prices <- object@pricePost
      mc     <- object@mcPost
      
      margins <- 1 - mc/prices
    }
    
    
    names(margins) <- object@labels
    
    return(as.vector(margins))
  }
  
)


#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "Auction2ndCap",
  definition=function(object,preMerger=TRUE,exAnte=TRUE){

    result <- calcProducerSurplus(object,preMerger=preMerger,exAnte=exAnte)/calcPrices(object,preMerger=preMerger,exAnte=exAnte)
    return(result)
  })

## compute margins
#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "Cournot",
  definition=function(object,preMerger=TRUE){


    if(preMerger){
      prices <- object@pricePre
    }
    else{prices <- object@pricePost}

    mc <- calcMC(object,preMerger = preMerger)
    prices <- matrix(prices, ncol=length(prices), nrow=length(mc),byrow=TRUE)



    margin <- 1 - mc/prices


    dimnames(margin) <- object@labels
    return(margin)
  }

)


#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "AIDS",
  definition=function(object,preMerger=TRUE){

    priceDelta <- object@priceDelta
    ownerPre   <- object@ownerPre
    shares     <- calcShares(object,TRUE)

    elastPre <-  t(elast(object,TRUE))
    marginPre <-  -1 * as.vector(MASS::ginv(elastPre * ownerPre) %*% (shares * diag(ownerPre))) / shares

    if(preMerger){
      names(marginPre) <- object@labels
      return(marginPre)}

    else{

      marginPost <- 1 - ((1 + object@mcDelta) * (1 - marginPre) / (priceDelta + 1) )
      names(marginPost) <- object@labels
      return(marginPost)
    }

  }
)

## compute margins
#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "LogitCap",
  definition=function(object,preMerger=TRUE){

    margins <- object@margins #capacity-constrained margins not identified -- use supplied margins

    if( preMerger) {
      capacities <- object@capacitiesPre

    }
    else{

      capacities <- object@capacitiesPost
    }



    quantities <- calcQuantities(object, preMerger=TRUE)
    constrained <-  abs(capacities - quantities) < 1e-5

    owner  <- object@ownerPre
    revenue<- calcShares(object,preMerger,revenue=TRUE)[!constrained]
    elast <-  elast(object,preMerger)
    margins[!constrained] <-  -1 * as.vector(MASS::ginv(t(elast*owner)[!constrained,!constrained]) %*% revenue) / revenue



    names(margins) <- object@labels

    return(as.vector(margins))
  }

)

## compute margins
#'@rdname Margins-Methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "Auction2ndLogit",
  definition=function(object,preMerger=TRUE,exAnte=FALSE,subset){


    nprods <- length(object@shares)

    if(missing(subset)){
      subset <- rep(TRUE,nprods)
    }

    if(!is.logical(subset) || length(subset) != nprods ){stop("'subset' must be a logical vector the same length as 'shares'")}

    margins <- rep(NA,nprods)

    if( preMerger) {
      owner  <- object@ownerPre}
    else{
      owner  <- object@ownerPost}


    owner <- owner[subset,subset]


    alpha <- object@slopes$alpha
    shares <- calcShares(object,preMerger=preMerger,revenue=FALSE)
    shares <- shares[subset]
    firmShares <- drop(owner %*% shares)
    margins[subset] <-  log(1-firmShares)/(alpha * firmShares)

    if(exAnte){ margins <-  margins * shares}

    names(margins) <- object@labels

    return(as.vector(margins))
  }

)

