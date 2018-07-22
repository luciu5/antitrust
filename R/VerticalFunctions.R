#' @title Supply Chain Merger Simulation
#' @name SupplyChain-Functions
#' @include BargainingClasses.R VerticalClasses.R
#' @aliases vertical
#' vertical.barg.bert
#' vertical.varg.2nd
#' @description Calibrates consumer demand using (Nested) Logit
#' and then simulates the price effect of a merger between two firms
#' under the assumption that all firms in the market are playing a
#' differentiated products Bertrand pricing game.
#' @description Let k denote the number of products produced by all
#' firms playing the Bertrand pricing game below.
#'
#' @param sharesDown A length k vector of product (quantity) shares. Values must be
#'   between 0 and 1.
#' @param pricesDown A length k vector of downstream product prices.
#' @param marginsDown A length k vector of downstream product margins, some of which may
#'   equal NA.
#' @param ownerPreDown EITHER a vector of length k whose values
#'   indicate which downstream firm produced a product pre-merger OR
#'   a k x k matrix of pre-merger ownership shares.
#' @param ownerPostDown EITHER a vector of length k whose values
#'   indicate which downstream firm produced a product after the merger OR
#'   a k x k matrix of post-merger ownership shares.
#' @param mcDeltaDown A vector of length k where each element equals the
#'   proportional change in a downstream firm's product-level marginal costs due to
#'     the merger. Default is 0, which assumes that the merger does not
#'     affect any products' marginal cost.
#' @param pricesUp A length k vector of upstream product prices.
#' @param marginsup A length k vector of upstream product margins, some of which may
#'   equal NA.
#' @param ownerPreUp EITHER a vector of length k whose values
#'   indicate which upstream firm produced a product pre-merger OR
#'   a k x k matrix of pre-merger ownership shares.
#' @param ownerPostUp EITHER a vector of length k whose values
#'   indicate which upstream firm produced a product after the merger OR
#'   a k x k matrix of post-merger ownership shares.
#' @param mcDeltaUp A vector of length k where each element equals the
#'   proportional change in a upstream firm's product-level marginal costs due to
#'     the merger. Default is 0, which assumes that the merger does not
#'     affect any products' marginal cost.
#' @param insideSize An integer equal to total pre-merger units sold.
#'   If shares sum to one, this also equals the size of the market.
#' @param normIndex An integer equalling the index (position) of the
#'   inside product whose mean valuation will be normalized to 1. Default
#'   is 1, unless \sQuote{shares} sum to less than 1, in which case the default is
#'   NA and an outside good is assumed to exist.
#' @param subset A vector of length k where each element equals TRUE if
#' the product indexed by that element should be included in the
#' post-merger simulation and FALSE if it should be excluded.Default is a
#' length k vector of TRUE.
#' @param priceOutside A length 1 vector indicating the price of the
#' outside good. Default is 0.
#' @param priceStartDown A length k vector of starting values used to solve for
#' downstream equilibrium prices. Default is the \sQuote{pricesDown} vector.
#' @param priceStartUp A length k vector of starting values used to solve for
#' upstream equilibrium price. Default is the \sQuote{pricesUp} vector.
#' @param isMax If TRUE, checks to see whether computed price equilibrium
#' locally maximizes firm profits and returns a warning if not. Default is FALSE.
#' @param parmsStart For \code{logit.cap.alm}, a length-2 vector of starting values
#' used to solve for the price coefficient and outside share (in that order). For
#' \code{logit.nets}, rhe first element should
#' always be the price coefficient and the remaining elements should be
#' the nesting parameters. Theory requires the nesting parameters to be
#' greater than the price coefficient. If missing then the random
#' draws with the appropriate restrictions are employed.
#' @param control.slopes A list of  \code{\link{optim}}
#' control parameters passed to the calibration routine
#' optimizer (typically the \code{calcSlopes} method).
#' @param control.equ A list of  \code{\link[BB]{BBsolve}} control parameters
#' passed to the non-linear equation solver (typically the \code{calcPrices} method).
#' @param constrain. Specify calibration strategy for estimating bargaining power parameter. 
#' "global" (default) assumes bargaining parameter is the same across all participants,"pair" assumes 
#' that all wholesaler/retailer pairs have a distinct parameter,"wholesaler" assumes that each wholesaler's
#' parameter is identical across negotiations, "retailer" assumes that each 
#' retailer's parameter is identical across negotiations. 
#' @param labels A k-length vector of labels. Default is "Prod#", where
#' \sQuote{#} is a number between 1 and the length of \sQuote{prices}.
#' @param ... Additional options to feed to the \code{\link[BB]{BBsolve}}
#' optimizer used to solve for equilibrium prices.
#'
#' @details Using product prices, quantity shares and all of the
#' product margins from at least one firm, \code{logit} is able to
#' recover the price coefficient and product mean valuations in a
#' Logit demand model. \code{logit} then uses these
#' calibrated parameters to simulate a merger between two firms.
#'
#'
#' @return \code{vertical.barg.bert} returns an instance of class
#' \code{\linkS4class{VertBargBertLogit}}.
#' @author Charles Taragin \email{charles.taragin@usdoj.gov}
#' @references 
#'
#' @examples 
#' 
#' ## Verical supply with 2 upstream firms,
#' ## 2 downstream firms, each offering 
#' ## a single product.
#' 
#' shareDown <- c( 0.25, 0.17, 0.25, 0.31)
#' marginDown <- c(0.28, 0.32, 0.33, 0.38)
#' priceDown <- c( 49.2, 42.2, 53.3, 46.3)
#' ownerPreDown <- paste0("D",rep(c(1,2),2))
#' marginUp <- c(0.27, 0.18, 0.25, 0.48)
#' priceUp <- c(47.4, 69.6, 85.5, 45.6)
#' ownerPreUp <- paste0("U",rep(c(1,2),2))
#'
#'
#' ## Simulate an upstream horizontal merger
#' ownerPostUp <- rep("U1",length(ownerPreUp))
#'
#'
#' simres_up <- vertical.barg.bert(sharesDown =shareDown,
#' pricesDown = priceDown,
#' marginsDown = marginDown,
#' ownerPreDown = ownerPreDown,
#' ownerPostDown = ownerPreDown,
#' pricesUp = priceUp,
#' marginsUp = marginUp,
#' ownerPreUp = ownerPreUp,
#' ownerPostUp = ownerPostUp)
#' 
#' 
#' ## simulate a horizontal 
#' @include VerticalClasses.R
NULL

#'@rdname SupplyChain-Functions
#'@export


vertical.barg.bert <- function(sharesDown,
                  pricesDown,marginsDown,
                  ownerPreDown,ownerPostDown,
                  mcDeltaDown=rep(0,length(pricesDown)),
                  pricesUp,marginsUp,
                  ownerPreUp,ownerPostUp,
                  mcDeltaUp=rep(0,length(pricesUp)),
                  normIndex=ifelse(isTRUE(all.equal(sum(shares),1,check.names=FALSE)),1, NA),
                  subset=rep(TRUE,length(prices)),
                  insideSize = NA_real_,
                  priceOutside = 0,
                  priceStartDown = pricesDown,
                  priceStartUp = pricesUp,
                  isMax=FALSE,
                  constrain = c("global","pair","wholesaler","retailer"),
                  control.slopes,
                  control.equ,
                  labels= paste("Prod1",1:length(pricesUp)),
                  ...
){


  constrain <-  match.arg(constrain)
  
  ## Create  containers to store relevant data
  
  up <- new("Bargaining",prices=pricesUp,
              shares = sharesDown,
              margins=marginsUp,
              ownerPre=ownerPreUp,
              ownerPost=ownerPostUp,
              mcDelta=mcDeltaUp,
              priceStart = priceStartUp,
              labels=labels)
  
  ## Convert downstream ownership vectors to ownership matrices
  up@ownerPre  <- ownerToMatrix(up,TRUE)
  up@ownerPost <- ownerToMatrix(up,FALSE)
  
  ## Calculate upstream marginal costs
  up@mcPre <-  calcMC(up,TRUE)
  up@mcPost <- calcMC(up,FALSE)
  
  down <- new("Logit",prices=pricesDown, shares=sharesDown,
                margins=marginsDown,
                normIndex=normIndex,
                ownerPre=ownerPreDown,
                ownerPost=ownerPostDown,
                insideSize=insideSize,
                mcDelta=mcDeltaDown,
                subset=subset,
                priceOutside=priceOutside,
                priceStart=priceStartDown,
                shareInside=ifelse(isTRUE(all.equal(sum(sharesDown),1,check.names=FALSE,tolerance=1e-3)),1,sum(sharesDown)),
                labels=labels)
  
  
  ## Convert downstream ownership vectors to ownership matrices
  down@ownerPre  <- ownerToMatrix(down,TRUE)
  down@ownerPost <- ownerToMatrix(down,FALSE)
  
  ## Calculate downstream demand slope coefficients
  down <- calcSlopes(down)
  
  
  ## Calculate downstream marginal costs
  down@mcPre <-  calcMC(down,TRUE)
  down@mcPost <- calcMC(down,FALSE)
  

  ##save upstream and downstream data
  result <- new("VertBargBertLogit",
                up = up,
                down = down
                )
  
  if(!missing(control.slopes)){
    result@control.slopes <- control.slopes
  }
  if(!missing(control.equ)){
    result@control.equ <- control.equ
  }


  ## compute bargaining paramters
  result <- calcSlopes(result, constrain = constrain)
 
  ## Solve Non-Linear System for Price Changes
  
  resultsPre  <- calcPrices(result,preMerger=TRUE)
  resultsPost <- calcPrices(result,preMerger=FALSE,subset=subset)

  result@down@pricePre <- resultsPre$down
  result@down@pricePost <- resultsPost$down
  
  result@up@pricePre <- resultsPre$up
  result@up@pricePost <- resultsPost$up
  
  return(result)

}
