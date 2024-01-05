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
#' @param levels If TRUE calculates CMCR in levels rather than as a percentage of pre-merger costs. Default is FALSE.
#' @param rel A length 1 character vector indicating whether CMCR should be calculated
#' relative to pre-merger cost (``cost'') or pre-merger price (``price''), 
#' Default is ``cost''. Ignored when levels is TRUE.
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
  definition=function(object, market=FALSE,levels=FALSE,rel=c("cost","price")){

    rel=match.arg(rel)
    
    isParty <- rowSums( abs(object@ownerPost - object@ownerPre) ) > 0

    output <- object@output
    outSign <- ifelse(output,1,-1) 
    ##Compute pre-merger margins
    marginPre  <- calcMargins(object,preMerger=TRUE)


    ##compute post-merger margins evaluated at pre-merger prices
    object@ownerPre <- object@ownerPost
    marginPost <- calcMargins(object,preMerger=TRUE)

    cmcr <-  (marginPost - marginPre)
    
    if(rel=="cost"){
     if(output){ cmcr <- cmcr/(1 - marginPre)}
      else{cmcr <- cmcr/(1 + marginPre)}
    }
    names(cmcr) <- object@labels

    cmcr <- cmcr[isParty]

    if(levels){
      
      if(rel=="cost") {cmcr <- cmcr*object@mcPre}
      else{cmcr <- cmcr*object@pricePre}
      
      cmcr[!isParty & is.na(cmcr)] <- 0}
    else{cmcr <- cmcr * 100}
    
    if(market){
      sharePost <- calcShares(object, preMerger=FALSE, revenue =FALSE)
      if(all(is.na(sharePost))) sharePost <- calcShares(object, preMerger=FALSE, revenue = TRUE)
      sharePost <- sharePost[isParty]
      sharePost <- sharePost/sum(sharePost)

      cmcr <- sum( cmcr * sharePost, na.rm = TRUE )
    }

    return(cmcr)
  }
)

#'@rdname CMCR-Methods
#'@export
setMethod(
  f= "cmcr",
  signature= "Cournot",
  definition=function(object,market=TRUE,levels=FALSE,rel=c("cost","price")){

    rel <- match.arg(rel)
    owner <- object@ownerPre
    isParty <- rowSums( abs(object@ownerPost - object@ownerPre) ) > 0

    if(market){
      shares <- calcShares(object,preMerger=TRUE,revenue=FALSE)
    shares[is.na(shares)] <- 0
    shares <- owner %*% shares
    if(any(isParty)) shares <- unique(shares[isParty,,drop=FALSE])
    if(nrow(shares) == 1 ){shares <- rbind(shares,shares)}
    mktElast <- elast(object, preMerger= TRUE,market=TRUE)

    denom <--mktElast * colSums(shares)
    if(rel=="cost" && !levels){ denom <-  denom -  colSums(shares^2)}
    cmcr <- 2 * apply(shares,2,prod) / denom

    if(levels){
      
      cmcr <- cmcr*object@pricePre
      
    }
    
   
    }
    
    else{ margins<- calcMargins(object,preMerger=TRUE,level=FALSE)
          if(any(isParty)) margins <- unique(margins[isParty,,drop=FALSE])
          cmcr <- apply(margins,2,rev)
          if(rel=="cost" && !levels) {cmcr <- cmcr/(1-margins)}
          
          if(levels){
            
            cmcr <- t(t(cmcr)*object@pricePre)
            
          }
          
    }
    
    
    
    
    if(!levels){cmcr <- cmcr * 100}
    
    return(cmcr)
  })


#'@rdname CMCR-Methods
#'@export
setMethod(
  f= "cmcr",
  signature= "AIDS",
  definition=function(object, market= FALSE,rel=c("cost","price")){

    rel <- match.arg(rel)
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

    cmcr <- marginPost - marginPre
    if(rel=="cost") cmcr <- cmcr/(1 - marginPre)
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
  definition=function(object,market=FALSE,levels=FALSE,rel=c("cost","price"),...){

    rel <- match.arg(rel)
    
    isParty <- rowSums(abs(object@ownerPost-object@ownerPre)) > 0
    
    partyStart <- object@mcDelta[isParty]
    
    sharesPost <- calcShares(object,
                             preMerger=FALSE,
                             revenue=FALSE)
    
    result <- object@mcDelta <- rep(0, length(isParty))  
    
    
    
    minD <- function(mc){
      
    
      object@mcDelta[isParty] <- mc  
      thispricedelta <- calcPriceDelta(object)
      
      return(thispricedelta[isParty])
    } 
    

    thisCMCR <- BB::BBsolve(partyStart,minD,quiet=TRUE,...) 
  
    if(thisCMCR$convergence != 0 ){stop("CMCR may not have successfully converged. Minimization Routine Reports:\n",thisCMCR$message)}
     
    result[isParty] <- -thisCMCR$par
    
    names(result) <- object@labels
    
    if(!levels){
      if(rel=="cost"){result <- result/object@mcPre}
      else{result <- result/object@pricePre}
      result <- result*100
    }
    
    if(market){
      sharesPost <- sharesPost[isParty]
      result <- sum(result[isParty]*sharesPost/sum(sharesPost),na.rm=TRUE)  
      
    }
    
    else{result <- result[isParty]}
    
    return(result)
    
  }
)
