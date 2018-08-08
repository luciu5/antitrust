#' @title Methods For Calculating Price Delta
#' @name PriceDelta-Methods
#' @docType methods

#' @aliases calcPriceDelta
#' calcPriceDelta,ANY-method
#' calcPriceDelta,Antitrust-method
#' calcPriceDelta,AIDS-method
#' calcPriceDelta,Auction2ndLogit-method
#' calcPriceDelta,Cournot-method
#' calcPriceDelta,VertBargBertLogit-method
#'
#' @description For Antitrust, the method computes equilibrium price changes
#' due to a merger assuming that firms are playing a
#' Nash-Bertrand or Cournot game. This is a wrapper method for computing
#' the difference between pre- and post-merger equilbrium prices.
#' @description For AIDS, the method computes equilibrium price changes
#' due to a merger assuming that firms are playing a
#' Nash-Bertrand or Cournot game and LA-AIDS. This method calls a non-linear
#' equations solver to find a sequence of price changes that satisfy
#' the Bertrand FOCs.
#' @param object An instance of one of the classes listed above.
#' @param levels If TRUE, report results in levels. If FALSE, report results in percents. Default is FALSE.
#' @param market If TRUE, calculates (post-merger) share-weighted average of metric. Default is FALSE.
#' @param isMax If TRUE, uses numerical derivatives to determine if
#' equilibrium price vector is a local maximum. Default is FALSE.
#' @param subset A vector of length k where each element equals TRUE if
#' the product indexed by that element should be included in the
#' post-merger simulation and FALSE if it should be excluded.Default is a
#' length k vector of TRUE.
#' @param exAnte If \sQuote{exAnte} equals TRUE then the
#' \emph{ex ante} expected result for each firm is produced, while FALSE produces the
#' expected result conditional on each firm winning the auction. Default is FALSE.
#' @param ... Additional values that may be used to change the default values of the non-linear
#' equation solver.
#'
#' @include ParamsMethods.R
#' @keywords methods
NULL

setGeneric (
  name= "calcPriceDelta",
  def=function(object,...){standardGeneric("calcPriceDelta")}
)

## Method to compute price changes
#'@rdname PriceDelta-Methods
#'@export
setMethod(
  f= "calcPriceDelta",
  signature= "Antitrust",
  definition=function(object, levels = FALSE, market = FALSE, ...  ){

    pricePre  <- object@pricePre
    pricePost <- object@pricePost

    if(levels){priceDelta <- pricePost - pricePre}
    else{priceDelta <- pricePost/pricePre - 1}
    #names(priceDelta) <- object@labels

    if(market){
      shares <- calcShares(object, ...)
      shares <- shares/sum(shares,na.rm=TRUE)
      priceDelta <- sum(priceDelta*shares,na.rm=TRUE)
    }

    return(priceDelta)

  }
)

#'@rdname PriceDelta-Methods
#'@export
setMethod(
  f= "calcPriceDelta",
  signature= "VertBargBertLogit",
  definition=function(object, levels = FALSE, market = FALSE, ...  ){
    
    up <- object@up
    down <- object@down
    
    downDelta <- calcPriceDelta(down,levels = levels, market = market, ...)
    upDelta   <-  calcPriceDelta(up,levels = levels, market = FALSE, ...)
    
    if(market){
      shares <- calcShares(down, ...)
      shares <- shares/sum(shares,na.rm=TRUE)
      upDelta <- sum(upDelta*shares,na.rm=TRUE)
    }
   
    priceDelta <- list(up = upDelta,
                       down= downDelta)
    
    return(priceDelta)
    
  }
)


#'@rdname PriceDelta-Methods
#'@export
setMethod(
  f= "calcPriceDelta",
  signature= "AIDS",
  definition=function(object,isMax=FALSE,levels=FALSE,subset,market=FALSE,...){

    if(market) return(sum(object@priceDelta * calcShares(object, preMerger = FALSE),na.rm=TRUE))

    ownerPost <- object@ownerPost

    nprods <- length(object@shares)
    if(missing(subset)){subset <- rep(TRUE,nprods)}

    if(!is.logical(subset) || length(subset) != nprods ){stop("'subset' must be a logical vector the same length as 'shares'")}


    ##Define system of FOC as a function of priceDelta
    FOC <- function(priceDelta){

      object@priceDelta <- exp(priceDelta)-1

      sharePost <-  calcShares(object,FALSE)
      elastPost <-  t(elast(object,FALSE))
      marginPost <- calcMargins(object,FALSE)


      thisFOC <- sharePost*diag(ownerPost) + as.vector((elastPost*ownerPost) %*% (sharePost*marginPost))
      thisFOC[!subset] <- sharePost[!subset]
      return(thisFOC)

    }




    ## Find price changes that set FOCs equal to 0
    minResult <- BBsolve(object@priceStart,FOC,quiet=TRUE,control=object@control.equ)


    if(minResult$convergence != 0){warning("'calcPrices' nonlinear solver may not have successfully converged. 'BBsolve' reports: '",minResult$message,"'")}

    if(isMax){

      hess <- genD(FOC,minResult$par) #compute the numerical approximation of the FOC hessian at optimium
      hess <- hess$D[,1:hess$p]


      if(any(eigen(hess)$values>0)){warning("Hessian of first-order conditions is not positive definite. Price vector may not maximize profits. Consider rerunning 'calcPrices' using different starting values")}
    }

    deltaPrice <- exp(minResult$par)-1
    names(deltaPrice) <- object@labels

    if(levels){deltaPrice <- calcPrices(object,FALSE) - calcPrices(object,TRUE)}

    if(market){
      sharePost <-  calcShares(object,FALSE,...)
      sharePost <- sharePost/sum(sharePost, na.rm=TRUE)

      deltaPrice <- sum(deltaPrice*sharePost,na.rm=TRUE)

    }



    return(deltaPrice)
  }
)


#'@rdname PriceDelta-Methods
#'@export
setMethod(
  f= "calcPriceDelta",
  signature= "Auction2ndLogit",
  definition=function(object,exAnte=FALSE,levels=TRUE){

    subset <- object@subset

    mcDelta <- object@mcDelta

    if(exAnte){
      sharesPost <- calcShares(object, preMerger=FALSE)
      mcDelta <- mcDelta*sharesPost
    }
    else{sharesPost <- rep(1,length(subset))}


    result <- calcMargins(object, preMerger=FALSE,exAnte=exAnte) + mcDelta -
      calcMargins(object, preMerger=TRUE,exAnte=exAnte)

    if(!levels){result <- result/calcPrices(object,preMerger = TRUE, exAnte = exAnte )}

    names(result) <- object@labels
    return(result)
  }
)
