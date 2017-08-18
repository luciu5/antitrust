\name{defineMarketTools-methods}
\docType{methods}
\alias{HypoMonTest}
\alias{HypoMonTest,ANY-method}
\alias{HypoMonTest,Bertrand-method}
\alias{calcPricesHypoMon}
\alias{calcPricesHypoMon,ANY-method}
\alias{calcPricesHypoMon,AIDS-method}
\alias{calcPricesHypoMon,Linear-method}
\alias{calcPricesHypoMon,LogLin-method}
\alias{calcPricesHypoMon,Logit-method}
\alias{calcPricesHypoMon,LogitCap-method}
\alias{calcPricesHypoMon,Auction2ndLogit-method}
\alias{calcPricesHypoMon,Cournot-method}
\alias{calcPriceDeltaHypoMon}
\alias{calcPriceDeltaHypoMon,ANY-method}
\alias{calcPriceDeltaHypoMon,AIDS-method}
\alias{calcPriceDeltaHypoMon,Bertrand-method}
\alias{calcPriceDeltaHypoMon,Cournot-method}
\alias{diversionHypoMon}
\alias{diversionHypoMon,ANY-method}
\alias{diversionHypoMon,AIDS-method}
\alias{diversionHypoMon,Bertrand-method}
\title{Methods For Implementing The Hypothetical Monopolist Test}

\description{
An Implementation of the Hypothetical Monopolist
Test described in the 2010 Horizontal Merger Guidelines.
}

\usage{
\S4method{HypoMonTest}{Bertrand}(object,prodIndex,ssnip=.05,...)
\S4method{calcPricesHypoMon}{ANY}(object,prodIndex)
\S4method{calcPriceDeltaHypoMon}{ANY}(object,prodIndex)
\S4method{diversionHypoMon}{Bertrand}(object,prodIndex,...)
\S4method{diversionHypoMon}{AIDS}(object)
}

\arguments{
  \item{}{Let k denote the number of products produced by all firms playing the
    Bertrand pricing game.}
  \item{object}{An instance of one of the classes listed above.}
  \item{prodIndex}{A vector of product indices that are to be placed under the control
    of the Hypothetical Monopolist. }
  \item{ssnip}{ A number between 0 and 1 that
    equals the threshold for a \dQuote{Small but Significant and
      Non-transitory Increase in Price} (SSNIP). Default is .05, or 5\%.}
    \item{...}{Pass options to the optimizer used to solve for
      equilibrium prices. }
  }

   \details{
     \code{HypoMonTest} is an implementation of the Hypothetical Monopolist Test
     on the products indexed by \sQuote{prodIndex} for a \sQuote{ssnip}. The
     Hypothetical Monopolist Test determines whether a profit-maximizing
     Hypothetical Monopolist who controls the products indexed by
     \sQuote{prodIndex} would increase the price of at least one of the merging
     parties' products in \sQuote{prodIndex} by a
     small, significant, and non-transitory amount (i.e. impose a SSNIP).

      \code{calcPriceDeltaHypoMon} calculates the price changes
      relative to (predicted) pre-merger prices that a
      Hypothetical Monopolist would impose on the products indexed by
      \sQuote{prodIndex},
      holding the prices of products not
      controlled by the Hypothetical Monopolist fixed at pre-merger
      levels. With the exception of \sQuote{AIDS}, the
      \code{calcPriceDeltaHypoMon} for all the classes listed above
       calls \code{calcPricesHypoMon} to compute price
      levels. \code{calcPriceDeltaHypoMon} is in turn called by
      \code{HypoMonTest}.

      \code{diversionHypoMon} calculates the matrix of revenue
    diversions between all products included in the merger simulation,
    \emph{irrespective} of whether or not they are also included in
    \sQuote{prodIndex}. This matrix is useful for diagnosing whether or not a
    product not included in \sQuote{prodIndex} may have a higher revenue diversion
    either to or from a product included in \sQuote{prodIndex}. Note that the \sQuote{AIDS}
    \code{diversionHypoMon} method does not contain the \sQuote{prodIndex}
    argument, as AIDS revenue diversions are only a function of demand parameters.
  }

    \value{
      \code{HypoMonTest} returns TRUE if a profit-maximizing
      Hypothetical Monopolist who controls the products indexed by
      \sQuote{prodIndex} would increase the price of at least one of the merging
      parties' products in \sQuote{prodIndex} by a \sQuote{ssnip}, and
      FALSE otherwise. \code{HypoMonTest} returns an error if \sQuote{prodIndex}
      does not contain at least one of the merging parties
      products.

      \code{calcPriceDeltaHypoMon} returns a vector of proportional price changes for
      all products placed under the control of the Hypothetical
      Monopolist (i.e. all products indexed by \sQuote{prodIndex}).\\
      \code{calcPricesHypoMon} is identical, but for price levels.

      \code{diversionHypoMon} returns a k x k matrix of diversions,
      where element i,j is the diversion from product i to product j.

    }

\references{
  U.S. Department of Justice and Federal Trade Commission,
  \emph{Horizontal Merger Guidelines}. Washington DC: U.S.
  Department of Justice, 2010.
  \url{http://www.justice.gov/atr/public/guidelines/hmg-2010.html}
  (accessed July 29, 2011).
  }

