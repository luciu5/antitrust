#' @title Supply Chain Merger Simulation
#' @name SupplyChain-Functions
#' @include VerticalClasses.R
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
#' @param supplyDown A length 1 character vector that specifies whether the downstream
#' game that firms are playing is a Nash-Bertrand Pricing game ("bertrand'') or a 2nd score 
#' auction ("2nd"). Default is "bertrand".
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
#' @param marginsUp A length k vector of upstream product margins, some of which may
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
#' @param control.slopes A list of  \code{\link{optim}}
#' control parameters passed to the calibration routine
#' optimizer (typically the \code{calcSlopes} method).
#' @param control.equ A list of  \code{\link[BB]{BBsolve}} control parameters
#' passed to the non-linear equation solver (typically the \code{calcPrices} method).
#' @param constrain Specify calibration strategy for estimating bargaining power parameter. 
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
#' @author Charles Taragin \email{ctaragin@ftc.gov}
#' @references 
#' Sheu, Gloria,  Taragin, Charles (2018). 
#' "Simulating Mergers in a Vertical Supply Chain with Bargaining," EAG Discussions Papers 201804, Department of Justice, Antitrust Division.
#' @examples 
#' \donttest{
#' ## Verical supply with 2 upstream firms,
#' ## 2 downstream firms, each offering 
#' ## a single product.
#' 
#' shareDown <- c( 0.1293482, 0.1422541, 0.4631014, 0.2152962)
#' marginDown <- c( 0.2067533, 0.2572215, 0.3082511, 0.3539681)
#' priceDown <- c( 63.08158, 50.70465, 95.82960, 83.45267)
#' ownerPreDown <- paste0("D",rep(c(1,2),each=2))
#' marginUp <- c(0.5810900, 0.5331135, 0.5810900, 0.5331135)
#' priceUp <- c( 40.11427, 27.73734, 40.11427, 27.73734)
#' ownerPreUp <- paste0("U",rep(c(1,2),2))
#' priceOutSide <- 10
#'
#'
#' ## Simulate an upstream horizontal merger
#' ownerPostUp <- rep("U1",length(ownerPreUp))
#'
#'
#' simres_up <- vertical.barg(sharesDown =shareDown,
#' pricesDown = priceDown,
#' marginsDown = marginDown,
#' ownerPreDown = ownerPreDown,
#' ownerPostDown = ownerPreDown,
#' pricesUp = priceUp,
#' marginsUp = marginUp,
#' ownerPreUp = ownerPreUp,
#' ownerPostUp = ownerPostUp,
#' priceOutside = priceOutSide)
#' 
#' 
#' ## simulate a horizontal 
#' }
#' @include VerticalClasses.R
NULL

#'@rdname SupplyChain-Functions
#@export


vertical.barg <- function(supplyDown = c("bertrand","2nd"), sharesDown,
                  pricesDown,marginsDown,
                  ownerPreDown,ownerPostDown,
                  mcDeltaDown=rep(0,length(pricesDown)),
                  pricesUp,marginsUp,
                  ownerPreUp,ownerPostUp,
                  mcDeltaUp=rep(0,length(pricesUp)),
                  normIndex=ifelse(isTRUE(all.equal(sum(sharesDown),1,check.names=FALSE)),1, NA),
                  subset=rep(TRUE,length(pricesDown)),
                  insideSize = NA_real_,
                  priceOutside = 0,
                  priceStartDown = pricesDown,
                  priceStartUp = pricesUp,
                  isMax=FALSE,
                  constrain = c("global","pair","wholesaler","retailer"),
                  control.slopes,
                  control.equ,
                  labels= paste0("Prod",1:length(pricesUp)),
                  ...
){

  stop("Work in progress")
  
  
  constrain <-  match.arg(constrain)
  
  if(supplyDown =="bertrand"){downClass = "Logit"}
  else{downClass == "Auction2ndLogit"}
  
  ## Create  containers to store relevant data
  
 
  
  up <- new("Bargaining",
              prices=pricesUp,
              shares = sharesDown,
              subset=subset,
              margins=marginsUp,
              ownerPre=ownerPreUp,
              ownerPost=ownerPostUp,
              mcDelta=mcDeltaUp,
              priceStart = priceStartUp,
              labels=labels,
              constrain =constrain)
  
  ## Convert downstream ownership vectors to ownership matrices
  up@ownerPre  <- ownerToMatrix(up,TRUE)
  up@ownerPost <- ownerToMatrix(up,FALSE)
  
  
 
  down <- new(downClass,prices=pricesDown, shares=sharesDown,
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

 
  ## Calculate marginal costs
  mcPre <- calcMC(result,TRUE)
  mcPost <- calcMC(result,FALSE)
  
  result@down@mcPre <-  mcPre$down
  result@down@mcPost <- mcPost$down
  
  result@up@mcPre <-  mcPre$up
  result@up@mcPost <- mcPost$up
  
  
   
  ## Solve Non-Linear System for Price Changes

  resultsPre  <- calcPrices(result,preMerger=TRUE)
  resultsPost <- calcPrices(result,preMerger=FALSE,subset=subset)

  result@down@pricePre <- resultsPre$down
  result@down@pricePost <- resultsPost$down
  
  result@up@pricePre <- resultsPre$up
  result@up@pricePost <- resultsPost$up
  
  return(result)

}


