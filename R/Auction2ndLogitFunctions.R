#' @title 2nd Score Procurement Auction Model with (Nested) Logit Demand
#' @name Auction2ndLogit-Functions
#' @aliases auction2nd.logit 
#' auction2nd.logit.nests 
#' auction2nd.logit.alm
#'
#' @description Calibrates consumer demand using (Nested) Logit and then
#' simulates the price effect of a merger between two firms
#' under the assumption that all firms in the market are playing a
#' differentiated products 2nd score auction game.
#' @description Let k denote the number of products produced by all firms playing the
#' auction game below.
#' @param prices A length k vector of product prices.
#' @param shares A length k vector of product (quantity) shares. Values must be
#' between 0 and 1.
#' @param margins A length k vector of product margins (in levels, not percents), some of which may
#' equal NA.
#' @param nests A length k factor of product nests.
#' @param diversions A k x k matrix of diversion ratios with diagonal
#' elements equal to -1. Default is missing, in which case diversion
#' according to share is assumed.
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
#' @param mktElast a negative value indicating market elasticity. Default is NA.
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
#' @param mcDeltaOutside A length 1 vector indicating the change in the marginal cost of the
#' outside good. Default is 0.
#' @param parmsStart For \code{auction2nd.logit.alm}, a length 2 vector of starting values used to solve for
#' price coefficient and the share of the outside good. The first element should
#' always be the price coefficient and the second should be
#' the outside good. For \code{auction2nd.logit.nests}, a length n+1 vector of starting values used to solve for
#' price coefficient and the nesting parameters. The first element should
#' always be the price coefficient and the remaining elements should be
#' the nesting parameters.
#' @param constraint if TRUE, then the nesting parameters for all
#' non-singleton nests are assumed equal. If FALSE, then each
#' non-singleton nest is permitted to have its own value. Default is
#' TRUE.
#' @param control.slopes A list of  \code{\link{optim}}  control parameters passed
#' to the calibration routine optimizer (typically the \code{calcSlopes} method).
#' @param labels A k-length vector of labels. Default is "Prod#", where
#' \sQuote{#} is a number between 1 and the length of \sQuote{prices}.
#'
#' @details Using product prices, quantity shares and all of the
#' product margins from at least one firm, \code{auction2nd.logit} is able to
#' recover the price coefficient and product mean valuations in a
#' Logit demand model. \code{auction2nd.logit} then uses these
#' calibrated parameters to simulate a merger between two firms, under the assumption that firms are participating in a 2nd score procurement auction.
#'
#'
#' \code{auction2nd.logit.nests} is identical to \code{auction2nd.logit} except that it assumes
#' that products can be grouped into nests. Additional margin information is needed to 
#' identify th nesting parameters. 
#' 
#' \code{auction2nd.logit.alm} is identical to \code{auction2nd.logit} except that it assumes
#' that an outside product exists and uses additional margin
#' information to estimate the share of the outside good.
#'
#' @return \code{auction2nd.logit} returns an instance of \code{\linkS4class{Auction2ndLogit}},
#' a child class of \code{\linkS4class{Logit}}. \code{auction2nd.logit.nests} returns an instance of \code{\linkS4class{Auction2ndLogitNests}}.
#' \code{auction2nd.logit} returns an instance of \code{\linkS4class{Auction2ndLogitALM}}.
#'
#' @seealso \code{\link{logit}},\code{\link{logit.nests}} for simulating mergers under a Nash-Bertrand pricing game with Logit demand 
#' @author Charles Taragin \email{ctaragin+antitrustr@gmail.com}
#' @references Miller, Nathan (2014). \dQuote{Modeling the effects of mergers in procurement}
#' \emph{International Journal of Industrial Organization} , \bold{37}, pp. 201-208.
#'
#' @examples
#' ## Calibration and simulation results from a merger between firms 2 and 3
#' ## of a 4-firm market
#' ## Source: Miller 2014 backup materials http://www.nathanhmiller.org/research
#'
#' share = c(0.29,0.40,0.28,0.03)
#'
#' price = c(35.53,  154, 84.08, 53.16)*1e3
#' cost  =  c(NA, 101, NA, NA)*1e3
#'
#' ownerPre <- ownerPost <- diag(length(share))
#'
#' #Suppose products 2 and 3 merge
#' ownerPost[2,3] <- ownerPost[3,2] <- 1
#'
#' margin = price - cost
#'
#' result.2nd <- auction2nd.logit(price,share,margin,
#'                            ownerPre=ownerPre,ownerPost=ownerPost,normIndex=2)
#'
#'
#' print(result.2nd)
#' summary(result.2nd,revenue=FALSE)
#'
#' ##re-run without any price information except Firm 2
#'
#' price <- rep(NA_real_, length(price))
#'
#' result.noprice <- auction2nd.logit(price,share,margin,
#'                                    ownerPre=ownerPre,ownerPost=ownerPost,normIndex=2)
#'
#' print(result.noprice)
#' summary(result.noprice,revenue=FALSE)
#'
#' ##changing the units of prices and margins can yield dramatically different results 
#'
#' price = c(35.53,  154, 84.08, 53.16)
#' cost  =  c(NA, 101, NA, NA)
#' margin <- price - cost
#'
#' result.units <- auction2nd.logit(price,share,margin,
#'                                    ownerPre=ownerPre,ownerPost=ownerPost,normIndex=2)
#'
#' print(result.units)
#' summary(result.units,revenue=FALSE)
#'
#' ## Get a detailed description of the 'Auction2ndLogit' class slots
#' showClass("Auction2ndLogit")
#'
#' ## Show all methods attached to the 'Auction2ndLogit' Class
#' showMethods(classes="Auction2ndLogit")
#'
#' @include Auction2ndCapFunctions.R
NULL

#'@rdname Auction2ndLogit-Functions
#'@export
auction2nd.logit <- function(prices,shares,margins,
                             ownerPre,ownerPost,
                             normIndex=ifelse(isTRUE(all.equal(sum(shares),1,check.names=FALSE)),1, NA),
                             mcDelta=rep(0,length(prices)),
                             subset=rep(TRUE,length(prices)),
                             insideSize = NA_real_,
                             mcDeltaOutside=0,
                             control.slopes,
                             labels=paste("Prod",1:length(prices),sep="")
){

  if(missing(prices)){prices <- rep(NA_integer_, length(shares))}
  ## Create Auction2ndLogit  container to store relevant data
  result <- new("Auction2ndLogit",prices=prices, shares=shares,
                margins=margins,
                normIndex=normIndex,
                ownerPre=ownerPre,
                ownerPost=ownerPost,
                insideSize = insideSize,
                mcDelta=mcDelta,
                subset=subset,
                priceOutside=mcDeltaOutside,
                shareInside=ifelse(isTRUE(all.equal(sum(shares),1,check.names=FALSE)),1,sum(shares)),
                priceStart=rep(0,length(shares)),
                labels=labels,
                cls = "Auction2ndLogit")

  if(!missing(control.slopes)){
    result@control.slopes <- control.slopes
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



#'@rdname Auction2ndLogit-Functions
#'@export
auction2nd.logit.nests <- function(prices,shares,margins, 
                                   nests, diversions,
                             ownerPre,ownerPost,
                             normIndex=ifelse(isTRUE(all.equal(sum(shares),1,check.names=FALSE)),1, NA),
                             mcDelta=rep(0,length(prices)),
                             subset=rep(TRUE,length(prices)),
                             insideSize = NA_real_,
                             mcDeltaOutside=0,
                             parmsStart,
                             constraint=TRUE,
                             control.slopes,
                             labels=paste("Prod",1:length(prices),sep="")
){
  
  if(missing(prices)){prices <- rep(NA_integer_, length(shares))}
  
  
  if(missing(diversions)){diversions <- matrix(NA,nrow=length(shares),ncol=length(shares))}
  
  
  nests <- factor(nests,levels = unique(nests)) # factor nests, keeping levels in the order in which they appear
  nNestParm <- sum(tapply(nests,nests,length)>1) # count the number of  non-singleton nests
  nMargins  <- length(margins[!is.na(margins)])
  maxNests  <- nMargins - 1
  
  
  if(nNestParm > maxNests){
    stop("Additional margins must be supplied in order to calibrate nesting parameters")
  }
  
  if(missing(parmsStart)){
    
    nNests <- nlevels(nests)
    parmsStart <- runif(nNests+1) # nesting parameter values are assumed to be between 0 and 1
    parmsStart[1] <- -1* parmsStart[1] # price coefficient is assumed to be negative
    
    if(constraint){parmsStart <- parmsStart[1:2]}
  }
  
  
  if(constraint && length(parmsStart)!=2){
    stop("when 'constraint' is TRUE, 'parmsStart' must be a vector of length 2")
  }
  else if(!constraint && nNestParm + 1 != length(parmsStart)){
    stop("when 'constraint' is FALSE, 'parmsStart' must be a vector of length ",nNestParm + 1)
    
  }
 
  ## Create Auction2ndLogitNests  container to store relevant data
  result <- new("Auction2ndLogitNests",prices=prices, shares=shares,
                margins=margins,
                normIndex=normIndex,
                nests=nests,
                diversion=diversions,
                ownerPre=ownerPre,
                ownerPost=ownerPost,
                insideSize = insideSize,
                mcDelta=mcDelta,
                subset=subset,
                priceOutside=mcDeltaOutside,
                shareInside=ifelse(isTRUE(all.equal(sum(shares),1,check.names=FALSE)),1,sum(shares)),
                priceStart=rep(0,length(shares)),
                parmsStart=parmsStart,
                constraint=constraint,
                labels=labels,
                cls = "Auction2ndLogit")
  
  if(!missing(control.slopes)){
    result@control.slopes <- control.slopes
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


#'@rdname Auction2ndLogit-Functions
#'@export
auction2nd.logit.alm <- function(prices,shares,margins,
                                 ownerPre,ownerPost,
                                 mktElast = NA_real_,
                                 insideSize = NA_real_,
                                 mcDelta=rep(0,length(prices)),
                                 subset=rep(TRUE,length(prices)),
                                 mcDeltaOutside=0,
                                 parmsStart,
                                 control.slopes,
                                 labels=paste("Prod",1:length(prices),sep="")
){


  if(missing(parmsStart)){
    parmsStart <- rep(.1,2)
    nm <- which(!is.na(margins))[1]
    parmsStart[1] <- 1/(margins[nm]*log(1-shares[nm])) #ballpark alpha for starting values
  }


  if(missing(prices)){prices <- rep(NA_integer_, length(shares))}


  ## Create Auction2ndLogitALM  container to store relevant data
  result <- new("Auction2ndLogitALM",prices=prices, shares=shares,
                margins=margins,
                ownerPre=ownerPre,
                ownerPost=ownerPost,
                mktElast = mktElast,
                insideSize = insideSize,
                mcDelta=mcDelta,
                subset=subset,
                priceOutside=mcDeltaOutside,
                shareInside=ifelse(isTRUE(all.equal(sum(shares),1,check.names=FALSE,tolerance=1e-3)),1,sum(shares)),
                priceStart=rep(0,length(shares)),
                parmsStart = parmsStart,
                labels=labels,
                cls = "Auction2ndLogit")

  if(!missing(control.slopes)){
    result@control.slopes <- control.slopes
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

