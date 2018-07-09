#' @title Methods For Calculating Compensating Marginal Cost Reductions
#' @name CMCR-Methods
#' @docType methods

#' @aliases cmcr-methods
#' cmcr,ANY-method
#' cmcr,Bertrand-method
#' cmcr,Auction2ndLogit-method
#' cmcr,Cournot-method
#' cmcr,AIDS-method
#'
#' @description   Calculate the marginal cost reductions necessary to restore
#' premerger prices in a merger, or the Upwards Pricing Pressure Index
#' for the products of merging firms playing a differentiated
#' products Bertrand pricing game.
#'
#' @param object An instance of one of the classes listed above.
#' @param market If TRUE, calculates (post-merger) share-weighted average of metric. Default is FALSE.
#' @param ... Additional arguments to pass to \code{cmcr}.
#'
#' @details \code{cmcr} uses the results from the merger simulation and calibration
#' methods associates with a particular class to compute the compensating
#' marginal cost reduction (CMCR) for each of the merging parties' products.
#'
#' @return \code{cmcr} returns a vector of length k equal to CMCR for the
#' merging parties' products and 0 for all other products.
#'
#' @seealso \code{\link{cmcr.bertrand}} is a function that calculates CMCR
#' without the need to first calibrate a demand system and simulate a
#' merger.
#'
#' @include ElastMethods.R
#' @keywords methods
NULL

setGeneric (
  name= "cmcr",
  def=function(object,...){standardGeneric("cmcr")}
)


#'@rdname CMCR-Methods
#'@export
##Method to compute Compensating Marginal Cost Reduction
setMethod(
  f= "cmcr",
  signature= "Bertrand",
  definition=function(object, market=FALSE){

    isParty <- rowSums( abs(object@ownerPost - object@ownerPre) ) > 0

    ##Compute pre-merger margins
    marginPre  <- calcMargins(object,preMerger=TRUE)


    ##compute post-merger margins evaluated at pre-merger prices
    object@ownerPre <- object@ownerPost
    marginPost <- calcMargins(object,preMerger=TRUE)

    cmcr <- (marginPost - marginPre)/(1 - marginPre)
    names(cmcr) <- object@labels

    cmcr <- cmcr[isParty]

    if(market){
      sharePost <- calcShares(object, preMerger=FALSE, revenue =FALSE)
      if(all(is.na(sharePost))) sharePost <- calcShares(object, preMerger=FALSE, revenue = TRUE)
      sharePost <- sharePost[isParty]
      sharePost <- sharePost/sum(sharePost)

      cmcr <- sum( cmcr * sharePost )
    }

    return(cmcr * 100)
  }
)

#'@rdname CMCR-Methods
#'@export
setMethod(
  f= "cmcr",
  signature= "Cournot",
  definition=function(object,...){

    owner <- object@ownerPre
    isParty <- rowSums( abs(object@ownerPost - object@ownerPre) ) > 0

    shares <- calcShares(object,preMerger=TRUE,revenue=FALSE)
    shares[is.na(shares)] <- 0
    shares <- owner %*% shares
    if(any(isParty)) shares <- unique(shares[isParty,,drop=FALSE])
    if(nrow(shares) == 1 ){shares <- rbind(shares,shares)}
    mktElast <- elast(object, preMerger= TRUE,market=TRUE)

    cmcr <- 2 * apply(shares,2,prod) / (-mktElast * colSums(shares) - colSums(shares^2) )

    return(cmcr * 100)
  })


#'@rdname CMCR-Methods
#'@export
setMethod(
  f= "cmcr",
  signature= "AIDS",
  definition=function(object, market= FALSE){


    ownerPre  <- object@ownerPre
    ownerPost <- object@ownerPost

    isParty <- rowSums( abs(ownerPost - ownerPre) ) > 0

    sharesPre <- calcShares(object,TRUE)
    sharesPre <- tcrossprod(1/sharesPre,sharesPre)

    marginPre <- calcMargins(object,TRUE)


    elastPre  <- t(elast(object,TRUE))

    divPre    <- elastPre/diag(elastPre)


    Bpost      <- divPre * sharesPre * ownerPost
    marginPost <- -1 * as.vector(MASS::ginv(Bpost) %*% (diag(ownerPost)/diag(elastPre))
    )

    cmcr <- (marginPost - marginPre)/(1 - marginPre)
    names(cmcr) <- object@labels

    cmcr <- cmcr[isParty]

    if(market){
      sharePost <- calcShares(object, preMerger=FALSE, revenue =FALSE)
      if(all(is.na(sharePost))) sharePost <- calcShares(object, preMerger=FALSE, revenue = TRUE)
      sharePost <- sharePost[isParty]
      sharePost <- sharePost/sum(sharePost)

      cmcr <- sum( cmcr * sharePost )
    }

    return(cmcr * 100)
  }
)


## CMCR does not make sense in 2nd score auction
## delete
#'@rdname CMCR-Methods
#'@export
setMethod(
  f= "cmcr",
  signature= "Auction2ndLogit",
  definition=function(object,...){

    stop("'cmcr' is currently not available")
  }
)
