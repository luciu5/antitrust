#' @title Auction Cap Methods
#' @name AuctionCap-Methods
#' @docType methods
#'
#' @aliases calcBuyerExpectedCost
#' calcBuyerValuation
#' calcExpectedLowestCost
#' calcExpectedPrice
#' calcOptimalReserve
#' calcSellerCostParms
#' cdfG
#' calcBuyerExpectedCost,Auction2ndCap-method
#' calcBuyerValuation,Auction2ndCap-method
#' calcExpectedLowestCost,Auction2ndCap-method
#' calcExpectedPrice,Auction2ndCap-method
#' calcOptimalReserve,Auction2ndCap-method
#' calcSellerCostParms,Auction2ndCap-method
#' cdfG,Auction2ndCap-method
#'
#' @param object An instance of the respective class (see description for the classes)
#' @param  preMerger If TRUE, the pre-merger ownership structure is used. If FALSE, the post-merger ownership structure is used.
#' Default is TRUE.
#' @param lower The minimum for the bidder's reserve price.
#' @param upper The maximum for the bidder's reserve price.
#' @param c \code{cdfG} calculates the probability that a cost draw less than or equal to
#' \sQuote{c} is realized for each firm. If \sQuote{c} is not supplied,
#' the buyer reserve and total capacity is used.
#' @param ... Additional arguments to pass to \code{calcSellerCostParms}
#'
#' @author Charles Taragin \email{ctaragin+antitrustr@gmail.com}
#' @examples showMethods(classes="Auction2ndCap") # show all methods defined for the class
#'
#' @description
#'   \code{calcBuyerExpectedCost} computes the expected amount that the buyer will pay
#'     to the auction winner.
#'
#'   \code{calcBuyerValuation} computes the value to the buyer of the outside option.
#'
#'   \code{calcExpectedLowestCost} computes the expected lowest cost of the winning bid.
#'
#'   \code{calcExpectedPrice} computes the expected price paid by the buyer.
#'
#'   \code{calcOptimalReserve} computes the bidder's optimal reserve price.
#'
#'   \code{calcSellerCostParms} calibrates the parameters of the Seller Cost
#'     CDF, as well as the reserve price, if not supplied.
#'
#'   \code{cdfG} calculates the probability that a cost draw less than or equal to
#'     \sQuote{c} is realized for each firm. If \sQuote{c} is not supplied,
#'     the buyer reserve and total capacity is used.
#'
#' @include ShowMethods.R
#' @keywords methods
NULL


setGeneric (
  name= "calcOptimalReserve",
  def=function(object,...){standardGeneric("calcOptimalReserve")}
)
setGeneric (
  name= "cdfG",
  def=function(object,...){standardGeneric("cdfG")}
)
setGeneric (
  name= "calcExpectedPrice",
  def=function(object,...){standardGeneric("calcExpectedPrice")}
)
setGeneric (
  name= "calcExpectedLowestCost",
  def=function(object,...){standardGeneric("calcExpectedLowestCost")}
)
setGeneric (
  name= "calcBuyerExpectedCost",
  def=function(object,...){standardGeneric("calcBuyerExpectedCost")}
)
setGeneric (
  name= "calcSellerCostParms",
  def=function(object,...){standardGeneric("calcSellerCostParms")}
)
setGeneric (
  name= "calcBuyerValuation",
  def=function(object,...){standardGeneric("calcBuyerValuation")}
)

#'@rdname AuctionCap-Methods
#'@export
setMethod(
  f= "calcSellerCostParms",
  signature= "Auction2ndCap",
  definition=function(object,...){

    ## method to calibrate seller cost distribution parameters
    sellerCostParms <- object@sellerCostParms
    reserve         <- object@reserve
    shareInside     <- object@shareInside
    margins         <- object@margins
    prices          <- object@prices
    parmsStart      <- object@parmsStart


    cdf <- object@sellerCostCDF

    minD <- function(parmsStart){

      if(is.na(reserve)){r <- parmsStart[1]; parmsStart <- parmsStart[-1]}
      else{r <- reserve}

      sellerCostParms <- parmsStart

      object@reservePre      <- r
      object@sellerCostParms <- sellerCostParms

      ## For uniform, frechet,  distribution bounds are a function
      ## of distribution parameters

      if(identical(cdf,"punif")){
        object@sellerCostBounds <- parmsStart }
      else if(identical(cdf,"pfrechet")){
        object@sellerCostBounds[1] <- parmsStart[1] }





      ##calculate each bidder's profit margin, conditional on bidder winning
      thisInShare <- cdfG(object,preMerger=TRUE)
      thisMargin  <- calcProducerSurplus(object,preMerger=TRUE,exAnte=FALSE)
      thisPrice   <- calcPrices(object,preMerger=TRUE,exAnte=FALSE)


      measure   <- c(margins - thisMargin/thisPrice,
                     1 - thisPrice/prices,
                     shareInside - thisInShare)

      return(sum(measure^2,na.rm=TRUE))

    }


    if(identical(cdf,"punif")) {

      ui = diag(length(parmsStart))
      if(is.na(reserve)){ui[1,nrow(ui)-1]=-1} #constrain reserve to be greater than cLower
      ui[nrow(ui),nrow(ui)-1]=-1 #constrain cLower to be less than cUpper
      ci = rep(0,length(parmsStart))


      result <- constrOptim(parmsStart,minD,grad=NULL,ui=ui,ci=ci,
                            control=object@control.slopes,...)

    }
    else{
      lb <- ub <- rep(Inf,length(parmsStart))

      if( identical(cdf,"pexp")){lb[1] <- 1e-20}
      else if( identical(cdf,"pweibull")){ lb[1:2] <- 1e-20} #shape,scale must be positive
      else if( identical(cdf,"pgumbel")){  lb[2] <- 1e-20} #scale must be positive
      else if( identical(cdf,"pfrechet")){ lb[2] <- 1e-20; lb[3] <- 2}  #scale must be positive, shape must be > 2 for finite variance

      if(is.na(reserve)){lb <- c(0,lb); ub <- c(Inf,ub)}

      if(length(parmsStart)>1){method="L-BFGS-B"}
      else{method="Brent"; ub=1e12} #'Brent' is equivalent to using optimize for 1D problems

      result <- optim(parmsStart,minD,method=method,lower=lb,upper=ub,
                      control=object@control.slopes,...)

    }

    result <- result$par
    if(is.na(reserve)){object@reserve <- result[1]; result <- result[-1]}
    object@sellerCostParms <- result

    if(identical(cdf, "punif")){
      object@sellerCostBounds <- result }
    else if(identical(cdf,"pfrechet")){
      object@sellerCostBounds[1] <- result[1] }


    return(object)
  }
)

#'@rdname AuctionCap-Methods
#'@export
setMethod(
  f= "calcBuyerValuation",
  signature= "Auction2ndCap",
  definition=function(object){

    ## Use FOC from buyers cost minimization problem
    ## to uncover buyer cost parameter
    capacities <- object@capacities
    totCap       <- sum(capacities)
    reserve    <- object@reserve

    object@reservePre <- reserve
    cdfF = match.fun(object@sellerCostCDF)
    pdfF = object@sellerCostPDF

    sellerCostParms <- c(list(reserve),as.list(object@sellerCostParms))
    fr = do.call(pdfF,sellerCostParms)

    sellerCostParms <- c(sellerCostParms,
                         lower.tail=as.list(object@sellerCostCDFLowerTail))
    Fr = do.call(cdfF,sellerCostParms)

    gr <- totCap*fr*(1-Fr)^(totCap-1)

    expectedPrice  <- calcExpectedPrice(object,preMerger=TRUE)
    partialSupplierProfits <- (1-Fr)^(totCap-capacities) - (1-Fr)^totCap
    partialSupplierProfits <- sum(partialSupplierProfits)/gr

    result <- reserve + partialSupplierProfits

    return(result)

  }
)

#'@rdname AuctionCap-Methods
#'@export
setMethod(
  f= "calcOptimalReserve",
  signature= "Auction2ndCap",
  definition=function(object,preMerger=TRUE,lower,upper){

    if(missing(lower)){lower <- max(object@sellerCostBounds[1],0)}
    if(missing(upper)){upper <- object@buyerValuation}

    minD <- function(r){
      if(preMerger){object@reservePre <- r}
      else{object@reservePost <- r}
      calcBuyerExpectedCost(object,preMerger=preMerger)
    }

    res <- optimize(
      f  = minD,
      lower = lower,
      upper = upper,
      tol=object@control.slopes$reltol
    )


    rStar <- res$minimum
    return(rStar)
  }
)

#'@rdname AuctionCap-Methods
#'@export
setMethod(
  f= "calcBuyerExpectedCost",
  signature= "Auction2ndCap",
  definition=function(object,preMerger=TRUE){


    shareInside <- cdfG(object,preMerger=preMerger)
    val  <- object@buyerValuation * (1-shareInside) + calcExpectedPrice(object,preMerger=preMerger)*shareInside

    return(val)
  })

#'@rdname AuctionCap-Methods
#'@export
setMethod(
  f= "cdfG",
  signature= "Auction2ndCap",
  definition=function(object,c,preMerger=TRUE){

    if(missing(c)){
      if(preMerger){
        c <- object@reservePre
        capacities <- sum(object@capacities)
      }
      else{
        c <- object@reservePost
        capacities <- sum(object@capacities*(1+object@mcDelta))
      }


    }

    else{

      if(preMerger){capacities <- object@capacities}
      else{capacities <- tapply(object@capacities*(1+object@mcDelta),object@ownerPost,sum)}


    }


    cdfF = match.fun(object@sellerCostCDF)
    sellerCostParms <- c(list(c),as.list(object@sellerCostParms),
                         lower.tail=as.list(object@sellerCostCDFLowerTail))

    Fc = do.call(cdfF,sellerCostParms)
    retval = 1-(1-Fc)^capacities

    if(!preMerger && length(capacities)>1){

      temp <- rep(NA, length(object@ownerPre))
      temp[object@ownerPre == object@ownerPost] <- retval
      retval <- temp

    }
    return(retval)

  }
)

#'@rdname AuctionCap-Methods
#'@export
setMethod(
  f= "calcExpectedPrice",
  signature= "Auction2ndCap",
  definition=function(object,preMerger=TRUE){

    val <- calcExpectedLowestCost(object,preMerger=preMerger) + sum(calcProducerSurplus(object,preMerger=preMerger),na.rm=TRUE)/cdfG(object,preMerger=preMerger)
    return(val)
  }
)

#'@rdname AuctionCap-Methods
#'@export
setMethod(
  f= "calcExpectedLowestCost",
  signature= "Auction2ndCap",
  definition=function(object,preMerger=TRUE){


    if(preMerger){capacities <- object@capacities}
    else{capacities <- object@capacities*(1+object@mcDelta)}

    num    <- calcMC(object,t=sum(capacities),preMerger=preMerger,exAnte=TRUE)


    retval <- num/cdfG(object, preMerger=preMerger)
    return(retval)


  })
