#' @title Compensating Marginal Cost Reductions and Upwards Pricing Pressure (Cournot)
#' @name CMCRCournot-Functions
#' @aliases cmcr.cournot 
#' cmcr.cournot2
#' upp.cournot
#' @description Calculate the marginal cost reductions necessary to restore
#' premerger prices (CMCR), or the net Upwards Pricing Pressure (UPP) in a
#' merger involving firms playing a 
#' homogeneous product Cournot pricing game.
#'
#' @param shares A length-2 vector containing merging party quantity shares.
#' @param mktElast A length-1 containing the industry elasticity.
#' @param prices A length-2 vector of product prices.
#' @param margins A length-2 vector of product margins.
#' @param party If TRUE calculate a length-2 vector of individial party CMCRs. 
#' If FALSE calculate share-weighted CMCR relative to share-weighted pre-merger
#' marginal costs. Default is FALSE
#' @param ownerPre EITHER a vector of length 2 whose values
#' indicate which of the merging parties produced a product pre-merger OR
#' a 2 x 2 matrix of pre-merger ownership shares.
#' @param ownerPost A 2 x 2 matrix of post-merger
#' ownership shares. Default is a 2 x 2 matrix of 1s.
#' @param mcDelta A vector of length 2 where each element equals the
#' proportional change in a product's marginal costs due to
#' the merger. Default is 0, which assumes that the merger does not
#' affect any products' marginal cost.
#' @param rel A length 1 character vector indicating whether CMCR should be calculated
#' relative to pre-merger cost (``cost'') or pre-merger price (``price''), 
#' Default is ``cost''.
#' @param labels A length-2 vector of product labels.
#'
#' @details The \sQuote{shares} (or \sQuote{margins}) vector must have 2 elements, 
#' and all \sQuote{shares} and \sQuote{margins}
#' elements must be between 0 and 1. The \sQuote{mktElast} vector must
#' have 1 non-negative  element.
#' @return when \sQuote{party} is FALSE (default), \code{cmcr.cournot}, \code{cmcr.cournot2} return a vector with 1 element whose value equals the percentage change
#' in the products' average marginal costs that the merged firms
#' must achieve in order to offset a price increase. When \sQuote{party} is TRUE,  \code{cmcr.cournot}, \code{cmcr.cournot2} return a vector with 2 element whose value equals the percentage change
#' in \emph{each} parties' marginal costs necessary to offset a price increase. When \sQuote{rel} equals "cost" (default) results are in terms of per-merger marginal costs. Otherwise, results are in terms of pre-merger price.
#' @references  Froeb, Luke and Werden, Gregory (1998).
#' \dQuote{A robust test for consumer welfare enhancing mergers among sellers
#' of a homogeneous product.}
#' \emph{Economics Letters}, \bold{58}(3), pp. 367 - 369.
#' 
#' Werden, Gregory and Froeb, Luke (2008). \dQuote{Unilateral Competitive Effects of Horizontal Mergers}, 
#' in Paolo Buccirossi (ed), Handbook of Antitrust Economics (MIT Press).
#' @author Charles Taragin
#' @seealso \code{\link{cmcr.bertrand}} for a differentiated products Bertrand version of this measure.
#'
#' @examples
#' shares=c(.05,.65)
#' industryElast = 1.9
#' margins=shares/industryElast
#'
#' ##  calculate average CMCR as a percentage of pre-merger costs
#' cmcr.cournot(shares,industryElast, rel="cost")
#' 
#' ##  calculate average CMCR as a percentage of pre-merger price
#' cmcr.cournot(shares,industryElast, rel="price")
#'
#' ##  calculate average CMCR using margins as a percentage of pre-merger costs
#' cmcr.cournot2(margins, party=TRUE,rel="cost")
#' 
#' ##  calculate the average CMCR for various shares and
#' ##  industry elasticities in a two-product merger where both firm
#' ##  products have identical share (see Froeb and
#' ##  Werden, 1998, pg. 369, Table 1)
#'
#'
#'
#' deltaHHI = c(100, 500, 1000, 2500, 5000) #start with change in HHI
#' shares = sqrt(deltaHHI/(2*100^2)) #recover shares from change in HHI
#' industryElast = 1:3
#'
#' result = matrix(nrow=length(deltaHHI),ncol=length(industryElast),
#'                 dimnames=list(deltaHHI,industryElast))
#'
#' for(s in 1:length(shares)){
#'   for(e in 1:length(industryElast)){
#'
#'
#'     result[s,e] = cmcr.cournot(rep(shares[s],2),industryElast[e])[1]
#'
#'   }}
#'
#' print(round(result,1))
#'
#' @include CMCRBertrandFunctions.R
NULL

#'@rdname CMCRCournot-Functions
#'@export
cmcr.cournot <- function(shares,mktElast,
                         party=FALSE,
                         rel=c("cost","price"),
                         labels=names(shares)){

  rel=match.arg(rel)
  
  if(!is.vector(shares) || length(shares)!=2) { stop("'shares'  must be a vector of length 2")}
  if(!is.vector(mktElast) || length(mktElast)!=1) { stop("'mktElast'  must be a vector of length 1")}
  if(any(shares < 0,na.rm=TRUE) ||
     any(shares >1,na.rm=TRUE) ) {
    stop("'shares'  must be between 0 and 1")}
  if(any(mktElast < 0,na.rm=TRUE)) { stop("'mktElast'  must be positive")}
  
  if(party){
  mcDelta <- rev(shares)/(mktElast-shares)  
  names(mcDelta) <- labels
  }
  
  else{
  denom <- mktElast*sum(shares)
  
  if(rel == "cost"){ denom <- denom - sum(shares^2)}
  
  mcDelta <- 2*prod(shares)/denom
}

  return(mcDelta*100)
}


#'@rdname CMCRCournot-Functions
#'@export
cmcr.cournot2 <- function(margins,
                          rel=c("cost","price"),
                          party=FALSE,
                          labels=names(margins)){
  
  rel=match.arg(rel)
  
  if(is.null(labels)){labels <- paste("Prod",1:length(margins),sep="")}
  if(!is.character(labels) || length(labels)!=length(margins)){
    stop("labels must be a  length 2 character vector")}
  
  if(!is.numeric(margins) || length(margins)!=2) { stop("'margins'  must be a vector of length 2")}
  if(any(margins < 0,na.rm=TRUE) ||
     any(margins >1,na.rm=TRUE) ) {
    stop("'margins'  must be between 0 and 1")}
  
  if(party){
  mcDelta <- rev(margins)
  
  if(rel =="cost"){ mcDelta <-  mcDelta/( 1- margins)}
  names(mcDelta) <- labels
  
  }
  else{
    mcDelta <- 2*prod(margins)
    
    denom <- sum(margins)
    
    if(rel =="cost"){ denom <-  denom - sum(margins^2)}
    
    mcDelta <- mcDelta/denom
  }
  
  
  
  return(mcDelta*100)
}

#'@rdname CMCRCournot-Functions
#'@export
upp.cournot <- function(prices, margins, ownerPre,
                        ownerPost=matrix(1,ncol=length(prices), nrow=length(prices)),
                        mcDelta=rep(0,length(prices)),
                        labels=names(margins))
{

  if(is.null(labels)){labels <- paste("Prod",1:length(prices),sep="")}
  
  if(!is.vector(prices) || length(prices) !=1){ stop("'prices'  must be a vector of length 1")}
  
  prices <- rep(prices,2)
  
  diversions <- rep(1,ncol=2,nrow=2); diag(diversions) <- -1

  result <- upp.bertrand(prices, margins, diversions, ownerPre,
                         ownerPost, mcDelta, labels)

  return(result)

}
