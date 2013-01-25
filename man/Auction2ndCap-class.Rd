\name{Auction2ndCap-class}
\Rdversion{1.0}
\docType{class}
\alias{Auction2ndCap-class}
%\alias{calcPrices,Auction2ndCap-method}
%\alias{calcShares,Auction2ndCap-method}
%\alias{calcSlopes,Auction2ndCap-method}
%\alias{CV,Auction2ndCap-method}
%\alias{elast,Auction2ndCap-method}

\title{Class "Auction2ndCap"}
\description{
   The "Auction2ndCap" class contains all the information needed to
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
    \item{\code{sellerCostCDF}:}{A function returning the Cumulative
      Distribution of Seller Costs}
    \item{\code{sellerCostPDF}:}{A function returning the Probability Density
      Distribution of Seller Costs}
    \item{\code{rangeCostDist}:}{The bounds of the seller's CDF and
      PDF}
    \item{\code{buyerCost}:}{Buyer's self-supply cost}
    \item{\code{buyerReserve}:}{Buyer's reservation price}
  }
}
\section{Extends}{
Class \code{\linkS4class{Antitrust}}, directly.
}

\section{Methods}{
  For all of methods containing the `preMerger' argument, `preMerger' takes
  on a value of TRUE or FALSE, where TRUE invokes the method using the
  pre-merger ownership structure, while FALSE invokes the method using
  the post-merger ownership structure.
  \describe{
    \item{\code{calcExpectedLowestCost}}{\code{signature(object = Auction2ndCap, preMerger
	= TRUE)}}{Computes the expected lowest cost of the winning bid.}
    \item{\code{calcExpectedSupplierProfits}}{\code{signature(object = Auction2ndCap, preMerger
	= TRUE)}}{}
    \item{\code{calcOptimalReserve}}{\code{signature(object = Auction2ndCap, preMerger
	= TRUE)}}{Computes Bidder's optimal reserve price.}
    \item{\code{calcPrices}}{\code{signature(object = Auction2ndCap, preMerger
	= TRUE)}}{Computes the expected price that the buyer pays,
      conditional on the buyer purchasing.}
    Hypothetical Monopolist would impose on its products
    relative to pre-merger prices. }
  \item{\code{calcBuyerExpectedCost}}{\code{signature(object = Auction2ndCap, preMerger
      = TRUE)}}{}
  \item{\code{cdfG}}{\code{signature(object = Auction2ndCap,
      c)}}{Calculates the probability that a cost draw less than or equal to
    `c' is realized. If `c' is not supplied, the buyer reserve is used.  }
}

\references{
%%  ~~put references to the literature/web site here~~
}
\author{Charles Taragin \email{charles.taragin@usdoj.gov}}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("Auction2ndCap")           # get a detailed description of the class
showMethods(classes="Auction2ndCap") # show all methods defined for the class
}
\keyword{classes}
