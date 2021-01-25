#' @title (Nested) Logit Demand Calibration and Merger Simulation)
#' @name Logit-Functions
#' @aliases logit
#' logit.nests
#' logit.nests.alm
#' logit.cap
#' logit.alm
#' logit.cap.alm
#' @description Calibrates consumer demand using (Nested) Logit
#' and then simulates the price effect of a merger between two firms
#' under the assumption that all firms in the market are playing a
#' differentiated products Bertrand pricing game.
#' @description Let k denote the number of products produced by all
#' firms playing the Bertrand pricing game below.
#'
#' @param prices A length k vector of product prices.
#' @param shares A length k vector of product (quantity) shares. Values must be
#'   between 0 and 1.
#' @param margins A length k vector of product margins, some of which may
#'   equal NA.
#' @param diversions A k x k matrix of diversion ratios with diagonal
#' elements equal to -1. Default is missing.
#' @param nests A length k vector identifying the nest that each
#'   product belongs to.
#' @param capacitiesPre A length k vector of pre-merger product capacities. Capacities
#'   must be at least as great as shares * insideSize.
#' @param capacitiesPost A length k vector of post-merger product capacities.
#' @param insideSize An integer equal to total pre-merger units sold.
#'   If shares sum to one, this also equals the size of the market.
#' @param normIndex An integer equalling the index (position) of the
#'   inside product whose mean valuation will be normalized to 1. Default
#'   is 1, unless \sQuote{shares} sum to less than 1, in which case the default is
#'   NA and an outside good is assumed to exist.
#' @param ownerPre EITHER a vector of length k whose values
#'   indicate which firm produced a product pre-merger OR
#'   a k x k matrix of pre-merger ownership shares.
#' @param ownerPost EITHER a vector of length k whose values
#'   indicate which firm produced a product after the merger OR
#'   a k x k matrix of post-merger ownership shares.
#' @param mktElast a negative value indicating market elasticity. Default is NA.
#' @param mcDelta A vector of length k where each element equals the
#'   proportional change in a product's marginal costs due to
#'     the merger. Default is 0, which assumes that the merger does not
#'     affect any products' marginal cost.
#' @param subset A vector of length k where each element equals TRUE if
#' the product indexed by that element should be included in the
#' post-merger simulation and FALSE if it should be excluded.Default is a
#' length k vector of TRUE.
#' @param constraint if TRUE, then the nesting parameters for all
#' non-singleton nests are assumed equal. If FALSE, then each
#' non-singleton nest is permitted to have its own value. Default is
#' TRUE.
#' @param priceOutside A length 1 vector indicating the price of the
#' outside good. Default is 0.
#' @param priceStart A length k vector of starting values used to solve for
#' equilibrium price. Default is the \sQuote{prices} vector.
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
#' \code{logit.alm} is identical to \code{logit} except that it assumes
#' that an outside product exists and uses additional margin
#' information to estimate the share of the outside good.
#' If market elasticity is known, it may be supplied using the
#' \sQuote{mktElast} argument.
#'
#' \code{logit.nests} is identical to \code{logit} except that it includes the \sQuote{nests}
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
#' \code{logit.nests.alm} is identical to \code{logit.nests} except that it assumes
#' that an outside product exists  and uses additional margin
#' information to estimate the share of the outside good.
#'
#' \code{logit.cap} is identical to  \code{logit}  except that firms are
#' playing the Bertrand pricing game under exogenously supplied capacity
#' constraints. Unlike \code{logit},  \code{logit.cap} requires users to
#' specify capacity constraints via \sQuote{capacities} and  the number of
#' potential customers in a market via \sQuote{mktSize}. \sQuote{mktSize} is needed to
#' transform \sQuote{shares} into quantities that must be directly compared to \sQuote{capacities}.
#'
#' In \code{logit}, \code{logit.nests} and \code{logit.cap},  if quantity shares sum to 1,
#' then one product's mean value is not identified and must be normalized
#' to 0. \sQuote{normIndex} may be used to specify the index (position) of the
#' product whose mean value is to be normalized. If the sum of revenue shares
#' is less than 1, both of these functions assume that the exists a k+1st
#' product in the market whose price and mean value are both normalized
#' to 0.
#'
#' @return \code{logit} returns an instance of class
#' \code{\linkS4class{Logit}}.
#' \code{logit.alm} returns an instance of \code{\linkS4class{LogitALM}}, a
#' child class of \code{\linkS4class{Logit}.}.
#' \code{logit.nests} returns an instance of \code{\linkS4class{LogitNests}}, a
#' child class of \code{\linkS4class{Logit}}.
#' \code{logit.cap} returns an instance of \code{\linkS4class{LogitCap}}, a
#' child class of \code{\linkS4class{Logit}.}
#'
#' @seealso \code{\link{ces}}
#' @author Charles Taragin \email{ctaragin@ftc.gov}
#' @references Anderson, Simon, Palma, Andre, and Francois Thisse (1992).
#' \emph{Discrete Choice Theory of Product Differentiation}.
#' The MIT Press, Cambridge, Mass.
#'
#'
#' Epstein, Roy and Rubinfeld, Daniel (2004).
#' \dQuote{Effects of Mergers Involving Differentiated Products.}
#'
#' Werden, Gregory and Froeb, Luke (1994).
#' \dQuote{The Effects of Mergers in
#'   Differentiated Products Industries: Structural Merger Policy and the
#'   Logit Model},
#' \emph{Journal of Law, Economics, \& Organization}, \bold{10}, pp. 407-426.
#'
#' Froeb, Luke, Tschantz, Steven  and Phillip Crooke (2003).
#' \dQuote{Bertrand Competition and Capacity Constraints: Mergers Among Parking Lots},
#' \emph{Journal of Econometrics}, \bold{113}, pp. 49-67.
#'
#' Froeb, Luke and Werden, Greg (1996).
#' \dQuote{Computational Economics and Finance: Modeling and Analysis with Mathematica, Volume 2.}
#' In Varian H (ed.), chapter Simulating Mergers among Noncooperative Oligopolists, pp. 177-95.
#' Springer-Verlag, New York.
#'
#' @examples ## Calibration and simulation results from a merger between Budweiser and
#' ## Old Style.
#' ## Source: Epstein/Rubenfeld 2004, pg 80
#'
#'
#' prodNames <- c("BUD","OLD STYLE","MILLER","MILLER-LITE","OTHER-LITE","OTHER-REG")
#' ownerPre <-c("BUD","OLD STYLE","MILLER","MILLER","OTHER-LITE","OTHER-REG")
#' ownerPost <-c("BUD","BUD","MILLER","MILLER","OTHER-LITE","OTHER-REG")
#' nests <- c("Reg","Reg","Reg","Light","Light","Reg")
#'
#' price    <- c(.0441,.0328,.0409,.0396,.0387,.0497)
#' shares   <- c(.066,.172,.253,.187,.099,.223)
#' margins <- c(.3830,.5515,.5421,.5557,.4453,.3769)
#'
#' insideSize <- 1000
#'
#' names(price) <-
#'   names(shares) <-
#'   names(margins) <-
#'   prodNames
#'
#' result.logit <- logit(price,shares,margins,
#'                       ownerPre=ownerPre,ownerPost=ownerPost,
#'                       insideSize = insideSize,
#'                       labels=prodNames)
#'
#'
#'
#' print(result.logit)           # return predicted price change
#' summary(result.logit)         # summarize merger simulation
#'
#' elast(result.logit,TRUE)      # returns premerger elasticities
#' elast(result.logit,FALSE)     # returns postmerger elasticities
#'
#' diversion(result.logit,TRUE)  # return premerger diversion ratios
#' diversion(result.logit,FALSE) # return postmerger diversion ratios
#'
#'
#' cmcr(result.logit)            #calculate compensating marginal cost reduction
#' upp(result.logit)            #calculate Upwards Pricing Pressure Index
#'
#' CV(result.logit)              #calculate representative agent compensating variation
#'
#'
#' ## Implement the Hypothetical Monopolist Test
#' ## for BUD and OLD STYLE using a 5\% SSNIP
#'
#' HypoMonTest(result.logit,prodIndex=1:2)
#'
#'
#'
#' ## Get a detailed description of the 'Logit' class slots
#' showClass("Logit")
#'
#' ## Show all methods attached to the 'Logit' Class
#' showMethods(classes="Logit")
#'
#' ## Show which classes have their own 'elast' method
#' showMethods("elast")
#'
#' ## Show the method definition for 'elast' and Class 'Logit'
#' getMethod("elast","Logit")
#'
#'
#'
#' #
#' # Logit With capacity Constraints
#' #
#'
#'
#' cap     <- c(66,200,300,200,99,300) # BUD and OTHER-LITE are capacity constrained
#' result.cap <- logit.cap(price,shares,margins,capacitiesPre=cap,
#'                         insideSize=insideSize,ownerPre=ownerPre,
#'                         ownerPost=ownerPost,labels=prodNames)
#' print(result.cap)
#'
#' @include LinearFunctions.R
NULL

#'@rdname Logit-Functions
#'@export
logit <- function(prices,shares,margins, diversions,
                  ownerPre,ownerPost,
                  normIndex=ifelse(isTRUE(all.equal(sum(shares),1,check.names=FALSE)),1, NA),
                  mcDelta=rep(0,length(prices)),
                  subset=rep(TRUE,length(prices)),
                  insideSize = NA_real_,
                  priceOutside = 0,
                  priceStart = prices,
                  isMax=FALSE,
                  control.slopes,
                  control.equ,
                  labels=paste("Prod",1:length(prices),sep=""),
                  ...
){

  if(missing(diversions)){diversions <- matrix(NA,nrow=length(shares),ncol=length(shares))}

  ## Create Logit  container to store relevant data
  result <- new("Logit",prices=prices, shares=shares,
                margins=margins, diversion = diversions,
                normIndex=normIndex,
                ownerPre=ownerPre,
                ownerPost=ownerPost,
                insideSize=insideSize,
                mcDelta=mcDelta,
                subset=subset,
                priceOutside=priceOutside,
                priceStart=priceStart,
                shareInside=ifelse(isTRUE(all.equal(sum(shares),1,check.names=FALSE,tolerance=1e-3)),1,sum(shares)),
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

#'@rdname Logit-Functions
#'@export
logit.nests <- function(prices,shares,margins, diversions,
                        ownerPre,ownerPost,
                        nests=rep(1,length(shares)),
                        normIndex=ifelse(sum(shares) < 1,NA,1),
                        mcDelta=rep(0,length(prices)),
                        subset=rep(TRUE,length(prices)),
                        priceOutside=0,
                        priceStart = prices,
                        isMax=FALSE,
                        constraint = TRUE,
                        parmsStart,
                        control.slopes,
                        control.equ,
                        labels=paste("Prod",1:length(prices),sep=""),
                        ...
){

  
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


  ## Create LogitNests  container to store relevant data
  result <- new("LogitNests",prices=prices, margins=margins, diversion = diversions,
                shares=shares,mcDelta=mcDelta,
                subset=subset,
                priceOutside=priceOutside,
                ownerPre=ownerPre,
                ownerPost=ownerPost,
                nests=nests,
                normIndex=normIndex,
                parmsStart=parmsStart,
                constraint=constraint,
                priceStart=priceStart,shareInside=sum(shares),
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

#'@rdname Logit-Functions
#'@export
logit.nests.alm <- function(prices,shares,margins,
                            ownerPre,ownerPost,
                            nests=rep(1,length(shares)),
                            mcDelta=rep(0,length(prices)),
                            subset=rep(TRUE,length(prices)),
                            priceOutside=0,
                            priceStart = prices,
                            isMax=FALSE,
                            constraint = TRUE,
                            parmsStart,
                            control.slopes,
                            control.equ,
                            labels=paste("Prod",1:length(prices),sep=""),
                            ...
){


  nests <- factor(nests,levels = unique(nests)) # factor nests, keeping levels in the order in which they appear
  nNestParm <- sum(tapply(nests,nests,length)>1) # count the number of  non-singleton nests
  nMargins  <- length(margins[!is.na(margins)])
  maxNests  <- nMargins - 2

  if(nNestParm > maxNests){
    stop("Additional margins must be supplied in order to calibrate nesting parameters")
  }

  if(missing(parmsStart)){

    nNests <- nlevels(nests)
    parmsStart <- runif(nNests+2) # nesting parameter values are assumed to be between 0 and 1
    parmsStart[1] <- -1* parmsStart[1] # price coefficient is assumed to be negative

    if(constraint){parmsStart <- parmsStart[1:3]}
  }


  if(constraint && length(parmsStart)!=3){
    stop("when 'constraint' is TRUE, 'parmsStart' must be a vector of length 3")
  }
  else if(!constraint && nNestParm + 2 != length(parmsStart)){
    stop("when 'constraint' is FALSE, 'parmsStart' must be a vector of length ",nNestParm + 2)

  }

  ## Create Logit  container to store relevant data
  result <- new("LogitNestsALM",prices=prices, shares=shares,
                margins=margins,
                ownerPre=ownerPre,
                ownerPost=ownerPost,
                nests=nests,
                mcDelta=mcDelta,
                subset=subset,
                priceOutside=priceOutside,
                priceStart=priceStart,
                shareInside=sum(shares),
                parmsStart=parmsStart,
                constraint=constraint,
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


#'@rdname Logit-Functions
#'@export
logit.cap <- function(prices,shares,margins,
                      ownerPre,ownerPost,
                      capacitiesPre = rep(Inf , length(prices)),
                      capacitiesPost = capacitiesPre,
                      insideSize,
                      normIndex=ifelse(sum(shares)<1,NA,1),
                      mcDelta=rep(0,length(prices)),
                      subset=rep(TRUE,length(prices)),
                      priceOutside=0,
                      priceStart = prices,
                      isMax=FALSE,
                      control.slopes,
                      control.equ,
                      labels=paste("Prod",1:length(prices),sep=""),
                      ...
){


  ## Create LogitCap  container to store relevant data
  result <- new("LogitCap",prices=prices, shares=shares,
                margins=margins,
                capacitiesPre=capacitiesPre,
                capacitiesPost=capacitiesPost,
                insideSize=insideSize,
                normIndex=normIndex,
                ownerPre=ownerPre,
                ownerPost=ownerPost,
                mcDelta=mcDelta,
                subset=subset,
                priceOutside=priceOutside,
                priceStart=priceStart,shareInside=sum(shares),
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


#'@rdname Logit-Functions
#'@export
logit.alm <- function(prices,shares,margins,
                      ownerPre,ownerPost,
                      mktElast = NA_real_,
                      insideSize = NA_real_,
                      mcDelta=rep(0,length(prices)),
                      subset=rep(TRUE,length(prices)),
                      priceOutside=0,
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
    parmsStart[1] <- -1/(margins[nm]*prices[nm]*(1-shares[nm])) #ballpark alpha for starting values
  }



  ## Create Logit  container to store relevant data
  result <- new("LogitALM",prices=prices, shares=shares,
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


#'@rdname Logit-Functions
#'@export
logit.cap.alm <- function(prices,shares,margins,
                          ownerPre,ownerPost,
                          capacitiesPre= rep(Inf,length(prices)),
                          capacitiesPost = capacitiesPre,
                          mktElast = NA_real_,
                          insideSize,
                          mcDelta=rep(0,length(prices)),
                          subset=rep(TRUE,length(prices)),
                          priceOutside=0,
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
    parmsStart[1] <- -1/(margins[nm]*prices[nm]*(1-shares[nm])) #ballpark alpha for starting values
  }


  ## Create LogitCap  container to store relevant data
  result <- new("LogitCapALM",prices=prices, shares=shares,
                margins=margins, 
                capacitiesPre=capacitiesPre,
                capacitiesPost=capacitiesPost,
                mktElast=mktElast,
                insideSize=insideSize,
                ownerPre=ownerPre,
                ownerPost=ownerPost,
                mcDelta=mcDelta,
                subset=subset,
                priceOutside=priceOutside,
                priceStart=priceStart,
                parmsStart=parmsStart,
                shareInside=ifelse(isTRUE(all.equal(sum(shares),1,check.names=FALSE,tolerance=1e-3)),1,sum(shares)),
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
