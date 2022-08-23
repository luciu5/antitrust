#' @title (Capacity Constrained) 2nd Price Auction Model
#' @name Auction2ndCap-Functions
#' @aliases Auction2ndCap
#' auction2nd.cap
#'
#' @description Calibrates the parameters of bidder cost distributions and
#' then simulates the price effect of a merger between two firms
#' under the assumption that firms are competing in a (Capacity
#' Constrained) 2nd price auction.
#' @description Let k denote the number of firms bidding in the auction below.
#'
#' @param capacities A length k vector of firm capacities OR capacity shares.
#' @param margins A length k vector of product margins. All margins must
#' be either be between 0 and 1, or NA.
#' @param prices A length k vector product prices. Prices may be NA.
#' @param reserve A length 1 vector equal to the buyer's reserve
#' price. Default is NA.
#' @param shareInside A length 1 vector equal to the probability that the
#' buyer does not select the outside option. Default is NA.
#' @param sellerCostCDF A length 1 character vector indicating which
#' probability distribution will be used  to model bidder cost
#' draws. Possible options are "punif", "pexp", "pweibull", "pgumbel",
#' "pfrechet". Default is "punif".
#' @param ownerPre A length k factor whose values
#' indicate which firms are present in the market pre-merger.
#' @param ownerPost A length k factor whose values
#' indicate which firms are present in the market post-merger.
#' @param mcDelta A vector of length k where each element equals the
#' proportional change in a firm's capacity  due to
#' the merger. Default is 0, which assumes that the merger does not
#' affect any products' capacity.
#' @param constrain.reserve If TRUE, the buyer's post-merger optimal
#' reserve price is assumed to equal the buyer's pre-merger optimal
#' reserve price. If FALSE, the buyer re-calculates her optimal reserve
#' price post-merger.
#' @param parmsStart A vector of starting values for calibrated parameters. See below
#' for more details.
#' @param control.slopes A list of  \code{\link{optim}}  control parameters passed to
#' the calibration routine optimizer (typically the \code{calcSlopes} method).
#' @param labels A k-length vector of labels. Default is "Firm", where
#' \sQuote{#} is a number between 1 and the length of \sQuote{capacities}.
#' @param ... Additional options to feed to either \code{\link{optim}} or \code{\link{constrOptim}}.
#'
#' @details \code{auction2nd.cap} examines how a merger affects equilbrium bidding
#' behavior when a single buyer is running a 2nd price procurement
#' auction with bidders whose marginal cost of supplying a homogenous product is
#' private information. This version of the model assumes that
#' bidders are differentiated by their capacities in the sense that firms
#' with greater capacity are more likely to have lower costs than firms
#' with smaller capacities.
#'
#' Using firm prices, shares, and margins, as well as information
#' on the auction reserve price as well as the proportion of buyers who
#' choose not to purchase from any bidder,  \code{auction2nd.cap}
#' calibrates the parameters of the common distribution from which
#' bidder's costs are drawn (and, if
#' not supplied, the implied reserve price) and then uses these calibrated
#' parameters to calibrate the value to the buyer of selecting the
#' outside option. Once these parameters have been calibrated,
#' \code{auction2nd.cap} computes the buyer's optimal pre-merger reservation
#' price, and if \sQuote{constrain.reserve} is FALSE, computes the buyer's
#' optimal post-merger reservation price (setting \sQuote{constrain.reserve} to
#' TRUE sets the buyer's post-merger optimal reserve equal to the buyer's
#' pre-merger optimal reserve). The pre- and post-merger expected price, conditional on a particular bidder winning, are then calculated.
#'
#' Currently, the common distribution from which costs may be drawn is
#' restricted to be either: Uniform ("punif"), Exponential ("pexp"), Weibull
#' ("pweibull"), Gumbel ("pgumbel"), or Frechet ("pfrechet"). Note that
#' the Exponential is a single parameter distribution, the Uniform and
#' Weibull are two parameter distributions, and the Gumbel and Frechet
#' are 3 parameter distributions. Accordingly, sufficient price, margin,
#' reserve, and outside share information must be supplied in order to
#' calibrate the parameters of the specified
#' distribution. \code{auction2nd.cap} returns an error if insufficient
#' information is supplied.
#'
#' @return \code{auction2nd.cap} returns an instance of class
#' \code{\linkS4class{Auction2ndCap}}.
#'
#' @author Charles Taragin \email{ctaragin+antitrustr@gmail.com}, with code
#' contributed by Michael Sandfort and Nathan Miller
#'
#' @references Keith Waehrer and  Perry, Martin (2003).
#' \dQuote{The Effects of Mergers in Open Auction Markets},
#' \emph{Rand Journal of Economics}, \bold{34(2)}, pp. 287-304.
#'
#' @examples
#' ##Suppose there are 3 firms (A,B,C) participating in a procurement auction with
#' ## an unknown reservation price and that firm A acquires firm B.
#'
#' caps <- c(0.65,0.30,0.05)           # total capacity normalized to 1 in this example
#' inShare    <- .67                   # probability that buyer does not select
#' # any bidder
#' prices     <- c(3.89, 3.79, 3.74)   # average price charged by each firm
#' margins    <- c(.228, .209, 0.197)  # average margin earned by each firm
#' ownerPre   <- ownerPost   <-c("A","B","C")
#' ownerPost[ownerPost=="B"] <- "A"
#'
#' ##assume costs are uniformly distributed with unknown bounds
#' result.unif = auction2nd.cap(
#'   capacities=caps,
#'   margins=margins,prices=prices,reserve=NA,
#'   shareInside=inShare,
#'   sellerCostCDF="punif",
#'   ownerPre=ownerPre,ownerPost=ownerPost,
#'   labels=ownerPre
#' )
#'
#' print(result.unif)
#' summary(result.unif)
#'
#' ## Get a detailed description of the 'Auction2ndCap' class slots
#' showClass("Auction2ndCap")
#'
#' ## Show all methods attached to the 'Auction2ndCap' Class
#' showMethods(classes="Auction2ndCap")
#'
#' @include AIDSFunctions.R
NULL


## Create Constructor Function
#'@rdname Auction2ndCap-Functions
#'@export
auction2nd.cap <- function(capacities, margins,prices,reserve=NA,shareInside=NA,
                           sellerCostCDF=c("punif","pexp","pweibull","pgumbel","pfrechet"),
                           ownerPre,ownerPost,
                           mcDelta=rep(0,length(capacities)),
                           constrain.reserve=TRUE, parmsStart,
                           control.slopes,
                           labels=as.character(ownerPre),...
){


  sellerCostCDF <- match.arg(sellerCostCDF)
  lower.tail    <- TRUE
  sellerCostPDF <- match.fun(paste("d",substring(sellerCostCDF,2),sep=""))

  if(is.na(reserve)){reserve <- NA_real_}
  if(is.na(shareInside)){shareInside <- NA_real_}

  if (missing(parmsStart)){
    avgPrice =  mean(prices,na.rm=TRUE)
    minPrice =  min(prices,na.rm=TRUE)
    maxPrice =  max(prices,na.rm=TRUE)


    if(identical(sellerCostCDF,"punif")){
      if(maxPrice > minPrice){
        parmsStart <- sellerCostBounds <- c(minPrice,maxPrice)} #uniform on price range
      else{parmsStart <- sellerCostBounds <- c(1e-20,maxPrice)}
    }
    else if(identical(sellerCostCDF,"pexp")){parmsStart=1/avgPrice; }
    else if(identical(sellerCostCDF,"pweibull")){ parmsStart=c(avgPrice,avgPrice)}
    else if(identical(sellerCostCDF,"pgumbel")){parmsStart=c(avgPrice,sqrt(avgPrice)) }
    else if(identical(sellerCostCDF,"pfrechet")){parmsStart=c(0,sqrt(avgPrice),2) }


    if(is.na(reserve) || missing(reserve)){
      parmsStart<-c(maxPrice,parmsStart)
    }
  }

  ##remove any names from parmStart
  names(parmsStart) <- NULL

  if(identical(sellerCostCDF,"punif")){sellerCostBounds <- parmsStart[-1]}
  else if(identical(sellerCostCDF,"pexp")){sellerCostBounds=c(0,Inf)}
  else if(identical(sellerCostCDF,"pweibull")){ sellerCostBounds=c(0,Inf)}
  else if(identical(sellerCostCDF,"pgumbel")){ sellerCostBounds=c(-Inf,Inf) ;lower.tail=FALSE}
  else if(identical(sellerCostCDF,"pfrechet")){ sellerCostBounds=c(parmsStart[2],Inf) ; lower.tail=FALSE}



  result <- new("Auction2ndCap",capacities=capacities,
                margins=margins,
                prices=prices,
                reserve=reserve,
                shareInside=shareInside,
                sellerCostCDF=sellerCostCDF,
                sellerCostCDFLowerTail=lower.tail,
                sellerCostBounds=sellerCostBounds,
                sellerCostPDF=sellerCostPDF,
                ownerPre=ownerPre,ownerPost=ownerPost,
                mcDelta=mcDelta,
                parmsStart=parmsStart,
                labels=labels)

  if(!missing(control.slopes)){
    result@control.slopes <- control.slopes
  }

  ## Calibrate seller cost parameters
  result                 <- calcSellerCostParms(result,...)

  ## Calibrate buyer cost parameter
  result@buyerValuation       <- calcBuyerValuation(result)

  ## Compute pre- and post-merger reserves
  result@reservePre      <- calcOptimalReserve(result,preMerger=TRUE) #Find Buyer Reserve pre-merger
  if(constrain.reserve){ result@reservePost <- result@reservePre}
  else{result@reservePost     <- calcOptimalReserve(result,preMerger=FALSE)} #Find Buyer Reserve post-merger

  ## Compute equilbrium prices
  result@pricePre         <- calcPrices(result,preMerger=TRUE,exAnte=FALSE)
  result@pricePost        <- calcPrices(result,preMerger=FALSE,exAnte=FALSE)

  ## Compute equilibrium marginal costs
  result@mcPre            <- calcMC(result,preMerger=TRUE,exAnte=FALSE)
  result@mcPost           <- calcMC(result,preMerger=FALSE,exAnte=FALSE)

  return(result)
}
