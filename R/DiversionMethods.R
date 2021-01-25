#' @title Methods For Calculating Diversion
#' @name Diversion-Methods
#' @docType methods
#'
#' @aliases diversion-methods
#' diversion
#' diversion,ANY-method
#' diversion,AIDS-method
#' diversion,Bertrand-method
#' diversion,VertBargBertLogit-method
#'
#' @description Calculate the diversion matrix between any two products in the market.
#' @param object An instance of one of the classes listed above.
#' @param preMerger If TRUE, calculates pre-merger price elasticities. If
#' FALSE, calculates post-merger price elasticities. Default is TRUE.
#' @param revenue If TRUE, calculates revenue diversion. If
#' FALSE, calculates quantity diversion. Default is TRUE for \sQuote{Bertrand}
#' and FALSE for \sQuote{AIDS}.
#' @return returns a k x k matrix of diversion ratios, where the i,jth
#' element is the diversion from i to j.
#'
#' @details For Bertrand, when \sQuote{revenue} is FALSE (the default),
#' this method uses the results from the merger calibration and
#' simulation to compute the \emph{quantity} diversion matrix between any two products
#' in the market. Element i,j of this matrix is the quantity diversion from
#' product i to product j, or the
#' proportion of product i's sales that leave (go to) i for (from) j due
#' to a increase (decrease) in i's price. Mathematically, quantity diversion is
#' \eqn{\frac{-\epsilon_{ji}share_j}{\epsilon_{ii}share_i}},
#' where \eqn{\epsilon_{ij}} is the cross-price elasticity from i to j.
#'
#' When \sQuote{revenue} is TRUE, this method computes the revenue diversion
#' matrix between any two products in the market. Element i,j of this matrix is the revenue diversion from
#' product i to product j, or the
#' proportion of product i's revenues that leave (go to) i for (from) j due
#' to a increase (decrease) in i's price. Mathematically, revenue diversion is
#' \eqn{-\frac{\epsilon_{ji}(\epsilon_{jj}-1)r_j}{\epsilon_{jj}(\epsilon_{ii}-1)r_j}}
#' where \eqn{r_i} is the revenue share of product i.
#'
#' When \sQuote{preMerger} is TRUE, diversions are
#' calculated at pre-merger equilibrium prices, and when \sQuote{preMerger} is FALSE, they
#' are calculated at post-merger equilibrium prices.
#'
#' @details For AIDS, when \sQuote{revenue} is TRUE (the default),
#' this method computes the \emph{revenue} diversion matrix between any two
#' products in the market. For AIDS, the revenue diversion from i to j is
#' \eqn{\frac{\beta_{ji}}{\beta_ij}}, where \eqn{\beta_{ij}} is the
#' percentage change in product i's revenue due to a change in j's price.
#'
#' When \sQuote{revenue} is FALSE, this \code{callNextMethod} is invoked. Will
#' yield a matrix of NAs if the user did not supply prices.
#'
#' When \sQuote{preMerger} is TRUE, diversions are
#' calculated at pre-merger equilibrium prices, and when \sQuote{preMerger} is FALSE, they
#' are calculated at post-merger equilibrium prices.
#'
#' @include UPPMethods.R
#' @keywords methods
NULL

setGeneric (
  name= "diversion",
  def=function(object,...){standardGeneric("diversion")}
)

## Method to compute diversion
#'@rdname Diversion-Methods
#'@export
setMethod(
  f= "diversion",
  signature= "Bertrand",
  definition=function(object,preMerger=TRUE,revenue=FALSE){



    labels  <- object@labels
    output  <-  calcShares(object,preMerger,revenue)


    elasticity <- elast(object,preMerger)



    if(revenue){
      diversion <- t(elasticity) / (diag(elasticity)-1)
      diversion <- -1 * diversion / diag(diversion)
    }

    else{
      diversion <- -1 * t(elasticity) / diag(elasticity)
    }

    diversion <- diversion * tcrossprod(1/output,output)
    dimnames(diversion) <-  list(labels,labels)

    return(diversion)
  }
)


#'@rdname Diversion-Methods
#'@export
setMethod(
  f= "diversion",
  signature= "AIDS",
  definition=function(object,preMerger=TRUE,revenue=TRUE){

    if(revenue){
      diversion <- -t(object@slopes)/diag(object@slopes)
      dimnames(diversion) <-  list(object@labels,object@labels)
      return(diversion)
    }

    else{callNextMethod(object,preMerger,revenue)}

  }

)


#'@rdname Diversion-Methods
#'@export
setMethod(
  f= "diversion",
  signature= "VertBargBertLogit",
  definition=function(object,preMerger=TRUE,revenue=TRUE){
    
    
    result <- diversion(object@down,preMerger=preMerger,revenue=revenue)
    
    return(result)
  })
