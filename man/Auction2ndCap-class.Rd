\name{Auction2ndCap-class}
\Rdversion{1.0}
\docType{class}
\alias{Auction2ndCap-class}
\alias{calcBuyerExpectedCost}
\alias{calcBuyerValuation}
\alias{calcExpectedLowestCost}
\alias{calcExpectedPrice}
\alias{calcOptimalReserve}
\alias{calcSellerCostParms}
\alias{cdfG}
\alias{calcBuyerExpectedCost,Auction2ndCap-method}
\alias{calcBuyerValuation,Auction2ndCap-method}
\alias{calcExpectedLowestCost,Auction2ndCap-method}
\alias{calcExpectedPrice,Auction2ndCap-method}
\alias{calcOptimalReserve,Auction2ndCap-method}
\alias{calcSellerCostParms,Auction2ndCap-method}
\alias{cdfG,Auction2ndCap-method}
\alias{calcMC,Auction2ndCap-method}
\alias{calcMargins,Auction2ndCap-method}
\alias{calcShares,Auction2ndCap-method}
\alias{calcPrices,Auction2ndCap-method}
\alias{calcProducerSurplus,Auction2ndCap-method}
\alias{summary,Auction2ndCap-method}



\title{Class \dQuote{Auction2ndCap}}
\description{
   The \dQuote{Auction2ndCap} class contains all the information needed to
   calibrate a 2nd price auction with capacity constraints
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by using the constructor function \code{\link{auction2nd.cap}}.
%%  ~~ describe objects here ~~
}

\section{Slots}{
  Let k denote the number of firms.
  \describe{
    \item{\code{capacities}:}{A length k vector of firm capacities.}
    \item{\code{margins}:}{A length k vector of product margins, some of which may
      equal NA.}
    \item{\code{prices}:}{A length k vector of product prices.}
    \item{\code{reserve}:}{A length 1 vector equal to observed buyer's
      reserve price. May equal NA.}
     \item{\code{shareInside}:}{A length 1 vector equal to the
       probability that a buyer does not select the outside option. May equal NA.}
    \item{\code{sellerCostCDF}:}{A length 1 character vector equal to
      the name of the function that calculates the Cumulative
      Distribution (CDF) of SellerCosts.}
    \item{\code{sellerCostCDFLowerTail}:}{A length 1 logical vector equal to
      TRUE if the  probabilities are P[X <= x], otherwise, P[X > x].}
    \item{\code{sellerCostPDF}:}{A function returning the Probability Density
      of Seller Costs.}
    \item{\code{sellerCostBounds}:}{The bounds on the seller's CDF.}
    \item{\code{sellerCostParms}:}{The parameters of the seller's CDF.}
    \item{\code{buyerValuation}:}{Buyer's self-supply cost.}
    \item{\code{reservePre}:}{Buyer's optimal pre-merger reservation price.}
    \item{\code{reservePost}:}{Buyer's optimal post-merger reservation
      price.}
    \item{\code{mcDelta}:}{A length k vector equal to the proportional
      change in a firm's capacity following the merger.}
    \item{\code{parmsStart}:}{A vector of starting values.}

  }
}
\section{Extends}{
Class \code{\linkS4class{Antitrust}}, directly.
}

\section{Methods}{
  For all of methods containing the \sQuote{preMerger} argument, \sQuote{preMerger} takes
  on a value of TRUE or FALSE, where TRUE invokes the method using the
  pre-merger ownership structure, while FALSE invokes the method using
  the post-merger ownership structure. Likewise, for all methods
  containing the \sQuote{exAnte} argument, if \sQuote{exAnte} equals TRUE then the
  \emph{ex ante} expected result for each firm is produced, while FALSE produces the
  expected result conditional on each firm winning the auction.
 \describe{
    \item{\code{calcBuyerExpectedCost}}{\code{signature(object = Auction2ndCap, preMerger = TRUE)}}{Computes the expected amount that the buyer will pay
  to the auction winner.}
    \item{\code{calcBuyerValuation}}{\code{signature(object = Auction2ndCap)}}{Computes the value to the buyer of the outside option.}
    \item{\code{calcExpectedLowestCost}}{\code{signature(object = Auction2ndCap, preMerger = TRUE)}}{Computes the expected lowest cost of the winning bid.}
    \item{\code{calcExpectedPrice}}{\code{signature(object = Auction2ndCap, preMerger = TRUE)}}{Computes the expected price paid by the buyer.}
    \item{\code{calcProducerSurplus}}{\code{signature(object = Auction2ndCap, preMerger = TRUE,exAnte=TRUE)}}{Computes the expected profits of each supplier}
    \item{\code{calcMC}}{\code{signature(object = Auction2ndCap, t, preMerger = TRUE,exAnte=TRUE)}}{Computes the expected marginal cost of
      each supplier for a given capacity profile \sQuote{t}. Default is \sQuote{preMerger} capacities.}
    \item{\code{calcMargins}}{\code{signature(object = Auction2ndCap, preMerger = TRUE,exAnte=TRUE)}}{Compute each firm's expected margin.}
    \item{\code{calcOptimalReserve}}{\code{signature(object = Auction2ndCap, preMerger = TRUE,lower,upper)}}{Computes the bidder's optimal reserve price.}
    \item{\code{calcPrices}}{\code{signature(object = Auction2ndCap, preMerger = TRUE,exAnte=TRUE)}}{Computes the expected price that the buyer pays,
      conditional on the buyer purchasing from a particular firm.}
    \item{\code{calcSellerCostParms}}{\code{signature(object = Auction2ndCap)}}{Calibrate the parameters of the Seller Cost
  CDF, as well as the reserve price, if not supplied.}
  \item{\code{calcShares}}{\code{signature(object = Auction2ndCap, preMerger = TRUE,exAnte=TRUE)}}{Compute the probability that a firm wins.}
  \item{\code{cdfG}}{\code{signature(object = Auction2ndCap,c,preMerger=TRUE)}}{Calculates the probability that a cost draw less than or equal to
    \sQuote{c} is realized for each firm. If \sQuote{c} is not supplied, the buyer reserve and total
    capacity is used.}
  \item{\code{summary}}{\code{signature(object =
      Auction2ndCap,exAnte=FALSE, parameters=FALSE,digits=2)}}{Summarize the
    results of the calibration and simulation.}
}}


\author{Charles Taragin \email{charles.taragin@usdoj.gov}}


\examples{
showClass("Auction2ndCap")           # get a detailed description of the class
showMethods(classes="Auction2ndCap") # show all methods defined for the class
}
\keyword{classes}
