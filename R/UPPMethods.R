#' @title Methods For Calculating Upwards Pricing Pressure Index (Bertrand)
#' @name UPP-Methods
#' @docType methods
#'
#' @aliases upp-methods
#' upp,ANY-method
#' upp,Bertrand-method
#' upp,AIDS-method
#' upp,Auction2ndLogit-method
#'
#' @description Calculate the Upwards Pricing Pressure Index
#' for the products of merging firms playing a differentiated
#' products Bertrand pricing game.
#'
#' @param object An instance of one of the classes listed above.
#'
#' @details \code{upp} uses the results from the merger simulation and
#' calibration to compute the upwards pricing pressure of the merger on
#' each merging parties' products.
#' @return \code{upp} returns a vector of length k equal to the  net UPP for the
#' merging parties' products and 0 for all other products.
#'
#' @seealso \code{\link{upp.bertrand}} calculates net UPP
#' without the need to first calibrate a demand system and simulate a merger.
#'
#' @include HHIMethods.R
#' @keywords methods
NULL

setGeneric (
  name= "upp",
  def=function(object,...){standardGeneric("upp")}
)


##Method to compute upp
#'@rdname UPP-Methods
#'@export
setMethod(
  f= "upp",
  signature= "Bertrand",
  definition=function(object){

    isParty     <- rowSums( abs(object@ownerPost  - object@ownerPre) ) > 0

    ownerPre    <- object@ownerPre
    ownerPost   <- object@ownerPost


    elastPre       <- elast(object,preMerger=TRUE)
    pricesPre       <- object@pricePre
    sharesPre      <- calcShares(object,preMerger=TRUE)

    mcPre       <- object@mcPre
    mcPost      <- object@mcPost


    marginsPre      <-  1 - mcPre/pricesPre
    marginsPost      <- 1 - mcPost/pricesPre

    focPre  <-  sharesPre*diag(ownerPre) +(t(elastPre)*ownerPre)  %*% (sharesPre*marginsPre)
    focPost <-  sharesPre*diag(ownerPost)+(t(elastPre)*ownerPost) %*% (sharesPre*marginsPost)


    result <- as.vector(focPost-focPre) #Generalized Pricing Pressure

    names(result) <- object@labels

    return(result[isParty])

  }
)

#'@rdname UPP-Methods
#'@export
setMethod(
  f= "upp",
  signature= "AIDS",
  definition=function(object){

    if(any(is.na(object@prices))){stop("UPP cannot be calculated without supplying values to 'prices'")}

    else{return(callNextMethod(object))}
  })


#'@rdname UPP-Methods
#'@export
setMethod(
  f= "upp",
  signature= "Auction2ndLogit",
  definition=function(object){

    shares <- calcShares(object,preMerger=TRUE,revenue=FALSE)
    diversion <- t(tcrossprod(shares,1/(1-shares)))
    diag(diversion) <- -1

    mcDelta <- object@mcDelta

    margins <- calcMargins(object,preMerger=TRUE)

    isParty <- abs(object@ownerPost - object@ownerPre)

    gross <- margins * shares * diversion * isParty
    gross <- sum(gross)

    isParty <- rowSums(isParty) > 0
    dpp <-  sum(mcDelta * diversion[,isParty], na.rm=TRUE)


    result <- gross + dpp


    return(result)

  }
)
