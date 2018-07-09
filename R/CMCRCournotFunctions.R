#' @title Compensating Marginal Cost Reductions and Upwards Pricing Pressure (Cournot)
#' @name CMCRCournot-Functions
#' @aliases cmcr.cournot
#' upp.cournot
#' @description Calculate the average marginal cost reduction necessary to restore
#' pre-merger prices, or the net Upwards Pricing Pressure in a
#' two-product merger involving firms playing a
#' homogeneous product Cournot pricing game.
#'
#' @param shares A length-2 vector containing merging party quantity shares.
#' @param mktElast A length-1 containing the industry elasticity.
#' @param prices A length-2 vector of product prices.
#' @param margins A length-2 vector of product margins.
#' @param ownerPre EITHER a vector of length 2 whose values
#' indicate which of the merging parties produced a product pre-merger OR
#' a 2 x 2 matrix of pre-merger ownership shares.
#' @param ownerPost A 2 x 2 matrix of post-merger
#' ownership shares. Default is a 2 x 2 matrix of 1s.
#' @param mcDelta A vector of length 2 where each element equals the
#' proportional change in a product's marginal costs due to
#' the merger. Default is 0, which assumes that the merger does not
#' affect any products' marginal cost.
#' @param labels A length-2 vector of product labels.
#'
#' @details The \sQuote{shares} vector must have 2 elements, and all \sQuote{shares}
#' elements must be between 0 and 1. The \sQuote{mktElast} vector must
#' have 1 non-negative  element.
#' @return A vector with 1 element whose value equals the percentage change
#' in the products' average marginal costs that the merged firms
#' must achieve in order to offset a price increase.
#' @references  Froeb, Luke and Werden, Gregory (1998).
#' \dQuote{A robust test for consumer welfare enhancing mergers among sellers
#' of a homogeneous product.}
#' \emph{Economics Letters}, \bold{58}(3), pp. 367 - 369.
#' @author Charles Taragin
#' @seealso \code{\link{cmcr.bertrand}} for a differentiated products Bertrand version of this measure.
#'
#' @examples
#' shares=c(.05,.65)
#' industryElast = 1.9
#'
#'
#' cmcr.cournot(shares,industryElast)
#'
#'
#' ## Calculate the necessary percentage cost reductions for various shares and
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
cmcr.cournot <- function(shares,mktElast){

  if(!is.vector(shares) || length(shares)!=2) { stop("'shares'  must be a vector of length 2")}
  if(!is.vector(mktElast) || length(mktElast)!=1) { stop("'mktElast'  must be a vector of length 1")}
  if(any(shares < 0,na.rm=TRUE) ||
     any(shares >1,na.rm=TRUE) ) {
    stop("'shares'  must be between 0 and 1")}
  if(any(mktElast < 0,na.rm=TRUE)) { stop("'mktElast'  must be positive")}

  mcDelta <- 2*prod(shares)/(mktElast*sum(shares) - sum(shares^2))


  return(mcDelta*100)
}

#'@rdname CMCRCournot-Functions
#'@export
upp.cournot <- function(prices, margins, ownerPre,
                        ownerPost=matrix(1,ncol=length(prices), nrow=length(prices)),
                        mcDelta=rep(0,length(prices)),
                        labels=paste("Prod",1:length(prices),sep=""))
{

  if(!is.vector(prices) || length(prices) !=2){ stop("'prices'  must be a vector of length 2")}
  diversions <- rep(1,ncol=2,nrow=2); diag(diversions) <- -1

  result <- upp.bertrand(prices, margins, diversions, ownerPre,
                         ownerPost, mcDelta, labels)

  return(result)

}
