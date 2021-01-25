#' @title Nash Bargaining Model with Logit Demand
#' @name BargainingLogit-Functions
#' @aliases bargaining.logit 
#' @description Calibrates consumer demand using  Logit and then
#' simulates the price effect of a merger between two firms
#' under the assumption that firms and customers in the market are playing a
#' differentiated products Nash Bargaining game.
#' @include Auction2ndLogitFunctions.R
#' @description Let k denote the number of products produced by all firms playing the
#' Nash Bargaining game below.
#' @param prices A length k vector of product prices.
#' @param shares A length k vector of product (quantity) shares. Values must be
#' between 0 and 1.
#' @param margins A length k vector of product margins (in levels, not percents), some of which may
#' equal NA.
# @param nests A length k factor of product nests.
# @param diversions A k x k matrix of diversion ratios with diagonal
# elements equal to -1. Default is missing, in which case diversion
# according to share is assumed.
#' @param normIndex An integer equalling the index (position) of the
#' inside product whose mean valuation will be normalized to 1. Default
#' is 1, unless \sQuote{shares} sum to less than 1, in which case the default is
#' NA and an outside good is assumed to exist.
#' @param ownerPre EITHER a vector of length k whose values
#' indicate which firm produced a product pre-merger OR
#' a k x k matrix of pre-merger ownership shares.
#' @param ownerPost EITHER a vector of length k whose values
#' indicate which firm produced a product after the merger OR
#' a k x k matrix of post-merger ownership shares.
#' @param bargpowerPre A length k vector of pre-merger bargaining power parameters. Values
#' must be between 0 (sellers have the power) and 1 (buyers the power). NA values are allowed,
#' though must be calibrated from additional margin and share data. Default is 0.5.
#' @param bargpowerPost A length k vector of post-merger bargaining power parameters. Values
#' must be between 0 (sellers have the power) and 1 (buyers the power). NA values are allowed,
#' though must be calibrated from additional margin and share data. Default is \sQuote{bargpowerPre}.
#' @param insideSize An integer equal to total pre-merger units sold.
#' If shares sum to one, this also equals the size of the market.
#' @param mcDelta A vector of length k where each element equals the
#' (level) change in a product's marginal costs due to
#' the merger. Default is 0, which assumes that the merger does not
#' affect any products' marginal cost.
#' @param subset A vector of length k where each element equals TRUE if
#' the product indexed by that element should be included in the
#' post-merger simulation and FALSE if it should be excluded.Default is a
#' length k vector of TRUE.
#' @param priceOutside A postive real number equal to the price of the outside good.
#' Default equals 1 for Logit demand.
#' @param priceStart A vector of length k who elements equal to an
#' initial guess of equilbrium prices. default is \sQuote{prices}.
#' @param control.slopes A list of  \code{\link{optim}}  control parameters passed
#' to the calibration routine optimizer (typically the \code{calcSlopes} method).
#' @param control.equ A list of  \code{\link[BB]{BBsolve}} control parameters passed
#' to the non-linear equation solver (typically the \code{calcPrices} method).
#' @param labels A k-length vector of labels. Default is "Prod#", where
#' \sQuote{#} is a number between 1 and the length of \sQuote{prices}.
#'
#' @details Using product prices, quantity shares and all of the
#' product margins from at least one firm, \code{auction2nd.logit} is able to
#' recover the price coefficient and product mean valuations in a
#' Logit demand model. \code{auction2nd.logit} then uses these
#' calibrated parameters to simulate a merger between two firms, under the assumption that firms are particpating in a 2nd score procurement auction.
#'
#'
#' @return \code{bargaining.logit} returns an instance of \code{\linkS4class{BargainingLogit}},
#' a child class of \code{\linkS4class{Logit}}.
#'
#' @seealso \code{\link{logit}} for simulating mergers under a Nash-Bertrand pricing game with Logit demand, and ,\code{\link{auction2nd.logit}}
#' for simulating mergers under a 2nd score auction with Logit demand. 
#' @author Charles Taragin \email{ctaragin@ftc.gov}
#' @references Miller, Nathan (2014). \dQuote{Modeling the effects of mergers in procurement}
#' \emph{International Journal of Industrial Organization} , \bold{37}, pp. 201-208.
#'
#' @examples
#' ## Calibration and simulation results from a merger between firms 2 and 3
#' ## of a 4-firm market
#' ## Source: Miller 2014 backup materials http://www.nathanhmiller.org/research
#'
#' share = c(0.29,0.40,0.28,0.03)
#' bargpower <- rep(0.6,4) # buyer has advantage
#' price = c(35.53,  154, 84.08, 53.16)
#' cost  =  c(NA, 101, NA, NA)
#'
#' ownerPre <- ownerPost <- diag(length(share))
#'
#' #Suppose products 2 and 3 merge
#' ownerPost[2,3] <- ownerPost[3,2] <- 1
#'
#' margin = (price - cost)/price
#'
#' result.barg <- bargaining.logit(price,share,margin,bargpowerPre=bargpower,
#'                            ownerPre=ownerPre,ownerPost=ownerPost,normIndex=2)
#'
#'
#' print(result.barg)
#' summary(result.barg,revenue=FALSE)
#'
#'
#' ## Get a detailed description of the 'BargainingLogit' class slots
#' showClass("BargainingLogit")
#'
#' ## Show all methods attached to the 'BargainingLogit' Class
#' showMethods(classes="BargainingLogit")
#'

#'@rdname BargainingLogit-Functions
#'@export
bargaining.logit <- function(prices,shares,margins,
                             ownerPre,ownerPost,
                             bargpowerPre=rep(0.5,length(prices)),
                             bargpowerPost=bargpowerPre,
                             normIndex=ifelse(isTRUE(all.equal(sum(shares),1,check.names=FALSE)),1, NA),
                             mcDelta=rep(0,length(prices)),
                             subset=rep(TRUE,length(prices)),
                             priceStart=prices,
                             insideSize = NA_real_,
                             priceOutside=0,
                             control.slopes,
                             control.equ,
                             labels=paste("Prod",1:length(prices),sep="")
){

  
  
  ## Create BargainingLogit  container to store relevant data
  result <- new("BargainingLogit",prices=prices, shares=shares,
                margins=margins,
                normIndex=normIndex,
                ownerPre=ownerPre,
                ownerPost=ownerPost,
                bargpowerPre=bargpowerPre,
                bargpowerPost=bargpowerPost,
                insideSize = insideSize,
                mcDelta=mcDelta,
                subset=subset,
                priceOutside=priceOutside,
                shareInside=ifelse(isTRUE(all.equal(sum(shares),1,check.names=FALSE)),1,sum(shares)),
                priceStart=priceStart,
                labels=labels,
                cls = "BargainingLogit")

  if(!missing(control.slopes)){
    result@control.slopes <- control.slopes
  }
  
  if(!missing(control.equ)){
    result@control.equ <- control.equ
  }
  
  ## Convert ownership vectors to ownership matrices
  result@ownerPre  <- ownerToMatrix(result,TRUE)
  result@ownerPost <- ownerToMatrix(result,FALSE)


  ## Calculate Demand Slope Coefficients
  result <- calcSlopes(result)

  ## Calculate marginal cost
  result@mcPre <-  calcMC(result,TRUE)
  result@mcPost <- calcMC(result,FALSE)

  
  ## Solve Non-Linear System for Price Changes
  result@pricePre  <- calcPrices(result,preMerger=TRUE)
  result@pricePost <- calcPrices(result,preMerger=FALSE)

  return(result)

}

