#' @title (Nested) Constant Elasticity of Substitution Demand Calibration and Merger Simulation
#' @name CES-Functions
#' @aliases ces
#' ces.alm
#' ces.nests
#'
#' @description Calibrates consumer demand using (Nested) Constant Elasticity of
#' Substitution (CES) and then simulates the price effect of a merger between two firms
#' under the assumption that all firms in the market are playing a
#' differentiated products Bertrand pricing game.
#' @description Let k denote the number of products produced by all firms playing the
#' Bertrand pricing game below.
#'
#' @param prices A length k vector of product prices.
#' @param shares A length k vector of product revenue shares.
#' @param margins A length k vector of product margins, some of which may
#' equal NA.
#' @param diversions A k x k matrix of diversion ratios with diagonal
#' elements equal to -1. Default is missing.
#' @param nests A length k vector identifying the nest that each
#' product belongs to.
#' @param ownerPre EITHER a vector of length k whose values
#' indicate which firm produced a product pre-merger OR
#' a k x k matrix of pre-merger ownership shares.
#' @param ownerPost EITHER a vector of length k whose values
#' indicate which firm produced a product after the merger OR
#' a k x k matrix of post-merger ownership shares.
#' @param mktElast a negative value indicating market elasticity. Default is NA.
#' @param normIndex An integer specifying the product index against which
#' the mean values of all other products are normalized. Default is 1.
#' @param insideSize total revenues included in the market.
#' @param mcDelta A vector of length k where each element equals the
#' proportional change in a product's marginal costs due to
#' the merger. Default is 0, which assumes that the merger does not
#' affect any products' marginal cost.
#' @param subset A vector of length k where each element equals TRUE if
#' the product indexed by that element should be included in the
#' post-merger simulation and FALSE if it should be excluded.Default is a
#' length k vector of TRUE.
#' @param constraint if TRUE, then the nesting parameters for all
#' non-singleton nests are assumed equal. If FALSE, then each
#' non-singleton nest is permitted to have its own value. Default is
#' TRUE.
#' @param priceOutside A length 1 vector indicating the price of the
#' outside good. Default is 1.
#' @param priceStart A length k vector of starting values used to solve for
#' equilibrium price. Default is the \sQuote{prices} vector.
#' @param isMax If TRUE, checks to see whether computed price equilibrium
#' locally maximizes firm profits and returns a warning if not. Default is FALSE.
#' @param parmsStart A vector of starting values used to solve for
#' price coefficient and nest parameters. The first element should
#' always be the price coefficient and the remaining elements should be
#' nesting parameters. Theory requires the nesting parameters to be
#' greater than the price coefficient. If missing then the random
#' draws with the appropriate restrictions are employed.
#' @param labels A k-length vector of labels. Default is "Prod#", where
#' \sQuote{#} is a number between 1 and the length of \sQuote{prices}.
#' @param control.slopes A list of  \code{\link{optim}}  control parameters passed to the
#' calibration routine optimizer (typically the \code{calcSlopes} method).
#' @param control.equ A list of  \code{\link[BB]{BBsolve}} control parameters passed
#' to the non-linear equation solver (typically the \code{calcPrices} method).
#' @param ... Additional options to feed to the \code{\link[BB]{BBsolve}}
#' optimizer used to solve for equilibrium prices.
#'
#' @details Using product prices, revenue shares and all of the
#' product margins from at least one firm, \code{ces} is able to
#' recover the price coefficient and product mean valuations in a
#' Constant Elasticity of Substitution demand model. \code{ces} then uses these
#' calibrated parameters to simulate the price effects of a merger between two firms under the
#' assumption that that all firms in the market are playing a
#' differentiated products Bertrand pricing game.
#'
#' \code{ces.alm} is identical to \code{ces} except that it assumes that an outside product exists and uses additional margin information to estimate the share of the outside good.
#'
#' \code{ces.nests} is identical to \code{ces} except that it includes the \sQuote{nests}
#' argument which may be used to assign products to different
#' nests. Nests are useful because they allow for richer substitution
#' patterns between products. Products within the same nest are assumed
#' to be closer substitutes than products in different nests. The degree
#' of substitutability between products located in different nests is
#' controlled by the value of the nesting parameter sigma.
#' The nesting parameters for singleton nests (nests containing
#'                                             only one product) are not identified and normalized to 1.  The vector of
#' sigmas is calibrated from the prices, revenue shares, and margins supplied
#' by the user.
#'
#' By default, all non-singleton nests are assumed to have a common value for sigma.
#' This constraint may be relaxed  by setting \sQuote{constraint} to
#' FALSE. In this case, at least one product margin must be supplied from
#' a product within each nest.
#'
#' In both  \code{ces} and  \code{ces.nests}, if revenue shares sum to 1,
#' then one product's mean value is not identified and must be normalized
#' to 1. \sQuote{normIndex} may be used to specify the index (position) of the
#' product whose mean value is to be normalized. If the sum of revenue shares
#' is less than 1, both of these functions assume that the exists a k+1st
#' product in the market whose price and mean value are both normalized
#' to 1.
#'
#' @return \code{ces} returns an instance of class \code{\linkS4class{CES}}.
#' \code{ces.alm} returns an instance of class \code{\linkS4class{CESALM}}.
#' \code{ces.nests} returns an instance of \code{\linkS4class{CESNests}}, a
#' child class of \code{\linkS4class{CES}}.
#' @seealso \code{\link{logit}}
#' @author Charles Taragin \email{ctaragin@ftc.gov}
#' @references Anderson, Simon, Palma, Andre, and Francois Thisse (1992).
#' \emph{Discrete Choice Theory of Product Differentiation}.
#' The MIT Press, Cambridge, Mass.
#'
#'
#' Epstein, Roy and Rubinfeld, Daniel (2004).
#' \dQuote{Effects of Mergers Involving Differentiated Products.}
#'
#' Sheu G (2011).
#' \dQuote{Price, Quality, and Variety: Measuring the Gains From Trade in Differentiated
#'   Products.}
#' U.S Department of Justice.
#'
#' @examples
#' ## Calibration and simulation results from a merger between Budweiser and
#' ## Old Style. Assume that typical consumer spends 1% of income on beer,
#' ## and that total beer expenditure in US is 1e9
#' ## Source: Epstein/Rubenfeld 2004, pg 80
#'
#' prodNames <- c("BUD","OLD STYLE","MILLER","MILLER-LITE","OTHER-LITE","OTHER-REG")
#' ownerPre <-c("BUD","OLD STYLE","MILLER","MILLER","OTHER-LITE","OTHER-REG")
#' ownerPost <-c("BUD","BUD","MILLER","MILLER","OTHER-LITE","OTHER-REG")
#' nests <- c("R","R","R","L","L","R")
#'
#' price    <- c(.0441,.0328,.0409,.0396,.0387,.0497)
#' shares   <- c(.071,.137,.251,.179,.093,.269)
#' margins  <- c(.3830,.5515,.5421,.5557,.4453,.3769)
#'
#' names(price) <-
#'   names(shares) <-
#'   names(margins) <-
#'   prodNames
#'
#' result.ces <-ces(price,shares,margins,ownerPre=ownerPre,ownerPost=ownerPost,
#'                  labels=prodNames)
#'
#' print(result.ces)           # return predicted price change
#' summary(result.ces)         # summarize merger simulation
#'
#' elast(result.ces,TRUE)      # returns premerger elasticities
#' elast(result.ces,FALSE)     # returns postmerger elasticities
#'
#' diversion(result.ces,TRUE)  # return premerger diversion ratios
#' diversion(result.ces,FALSE) # return postmerger diversion ratios
#'
#' cmcr(result.ces)            #calculate compensating marginal cost reduction
#' upp(result.ces)             #calculate Upwards Pricing Pressure Index
#'
#' CV(result.ces)              #calculate compensating variation as a percent of
#' #representative consumer income
#'
#' ## Implement the Hypothetical Monopolist Test
#' ## for BUD and OLD STYLE using a 5\% SSNIP
#'
#' HypoMonTest(result.ces,prodIndex=1:2)
#'
#'
#' ## Get a detailed description of the 'CES' class slots
#' showClass("CES")
#'
#' ## Show all methods attached to the 'CES' Class
#' showMethods(classes="CES")
#'
#' ## Show which class have their own 'elast' method
#' showMethods("elast")
#'
#' ## Show the method definition for 'elast' and Class 'CES'
#' getMethod("elast","CES")
#'
#' @include BertrandFunctions.R
NULL

#'@rdname CES-Functions
#'@export
ces <- function(prices,shares,margins,  diversions,
                ownerPre,ownerPost,
                normIndex=ifelse(sum(shares)<1,NA,1),
                insideSize = NA_real_,
                mcDelta=rep(0,length(prices)),
                subset=rep(TRUE,length(prices)),
                priceOutside=1,
                priceStart = prices,
                isMax=FALSE,
                control.slopes,
                control.equ,
                labels=paste("Prod",1:length(prices),sep=""),
                ...
){


  if(missing(diversions)){diversions <- matrix(NA,nrow=length(shares),ncol=length(shares))}

  ## Create CES  container to store relevant data
  result <- new("CES",prices=prices, shares=shares, margins=margins, 
                diversion=diversions,
                normIndex=normIndex,
                mcDelta=mcDelta,
                insideSize = insideSize,
                subset=subset,
                priceOutside=priceOutside,
                ownerPre=ownerPre,
                ownerPost=ownerPost,
                priceStart=priceStart,
                shareInside= sum(shares),labels=labels)

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
  result@pricePre  <- calcPrices(result,preMerger=TRUE,isMax=isMax,...)
  result@pricePost <- calcPrices(result,preMerger=FALSE,isMax=isMax,subset=subset,...)

  return(result)

}


#'@rdname CES-Functions
#'@export
ces.alm <- function(prices,shares,margins,
                    ownerPre,ownerPost,
                    mktElast = NA_real_,
                    insideSize = NA_real_,
                    mcDelta=rep(0,length(prices)),
                    subset=rep(TRUE,length(prices)),
                    priceOutside=1,
                    priceStart = prices,
                    isMax=FALSE,
                    parmsStart,
                    control.slopes,
                    control.equ,
                    labels=paste("Prod",1:length(prices),sep=""),
                    ...
){


  if(missing(parmsStart)){
    parmsStart <- rep(.1,2)
    nm <- which(!is.na(margins))[1]
    parmsStart[1] <- 1/(margins[nm]*(1-shares[nm])) - shares[nm]/(1-shares[nm]) #ballpark gamma for starting values
  }



  ## Create CES  container to store relevant data
  result <- new("CESALM",prices=prices, shares=shares,
                margins=margins,
                ownerPre=ownerPre,
                ownerPost=ownerPost,
                mktElast = mktElast,
                insideSize = insideSize,
                mcDelta=mcDelta,
                subset=subset,
                priceOutside=priceOutside,
                priceStart=priceStart,
                shareInside=ifelse(isTRUE(all.equal(sum(shares),1,check.names=FALSE,tolerance=1e-3)),1,sum(shares)),
                parmsStart=parmsStart,
                labels=labels)



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
  result@pricePre  <- calcPrices(result,preMerger=TRUE,isMax=isMax,...)
  result@pricePost <- calcPrices(result,preMerger=FALSE,isMax=isMax,subset=subset,...)

  return(result)

}


#'@rdname CES-Functions
#'@export
ces.nests <- function(prices,shares,margins,
                      ownerPre,ownerPost,
                      nests=rep(1,length(shares)),
                      normIndex=ifelse(sum(shares) < 1,NA,1),
                      insideSize = NA_real_,
                      mcDelta=rep(0,length(prices)),
                      subset=rep(TRUE,length(prices)),
                      priceOutside=1,
                      priceStart = prices,
                      isMax=FALSE,
                      constraint = TRUE,
                      parmsStart,
                      control.slopes,
                      control.equ,
                      labels=paste("Prod",1:length(prices),sep=""),
                      ...
){



  nests <- factor(nests,levels=unique(nests))


  if(missing(parmsStart)){
    nNests <- nlevels(nests)
    parmsStart <- cumsum(runif(nNests+1,1,1.5)) # parameter values are assumed to be greater than 1

    if(constraint){parmsStart <- parmsStart[1:2]}
  }

  ## Create CESNests  container to store relevant data
  result <- new("CESNests",prices=prices, shares=shares,margins=margins,
                mcDelta=mcDelta,
                insideSize = insideSize,
                subset=subset,
                priceOutside=priceOutside,
                ownerPre=ownerPre,
                ownerPost=ownerPost,
                nests=nests,
                normIndex=normIndex,
                parmsStart=parmsStart,
                priceStart=priceStart,
                constraint=constraint,
                shareInside=sum(shares),labels=labels)


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
  result@pricePre  <- calcPrices(result,preMerger=TRUE,isMax=isMax,...)
  result@pricePost <- calcPrices(result,preMerger=FALSE,isMax=isMax,subset=subset,...)


  return(result)

}




