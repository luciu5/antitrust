\name{defineMarkets-methods}
\docType{methods}
\alias{isMarket-methods}
\alias{smallestMarket-methods}
\alias{calcPriceDeltaHypoMon-methods}
\alias{isMarket}
\alias{smallestMarket}
\alias{isMarket,Bertrand-methods}
\alias{smallestMarket,Bertrand-methods}
\alias{calcPriceDeltaHypoMon}
\alias{calcPriceDeltaHypoMon,Linear-method}
\alias{calcPriceDeltaHypoMon,LogLin-method}
\alias{calcPriceDeltaHypoMon,Logit-method}
\alias{calcPriceDeltaHypoMon,AIDS-method}
\title{Methods For Defining Antitrust Markets}

\description{
Identify candidate antitrust markets using the Hypothetical Monopolist
Test described in the 2010 Horizontal Merger Guidelines.
}

\usage{
\S4method{isMarket}{Bertrand}(object,prodIndex,ssnip=.05,supplyResponse=FALSE,...)
\S4method{smallestMarket}{Bertrand}(object,startProd,prodOrder,ssnip=.05,supplyResponse=FALSE,...)

\S4method{calcPriceDeltaHypoMon}{Logit,AIDS,Linear,LogLin}(object,prodIndex)
}

\arguments{
  \item{}{Let k denote the number of products produced by all firms playing the
    Bertrand pricing game.}
  \item{object}{An instance of one of the classes listed above.}
  \item{prodIndex}{A vector of product indices that are to be placed under the control
    of the Hypothetical Monopolist. }
  \item{ssnip}{ A number, typically between 0 and 1 that
    equals the threshold for a "Small but Significant and
    Non-transitory Increase in Price" (SSNIP). Default is .05, or 5\%.}
  \item{supplyResponse}{ When FALSE (the default), sets the prices of
    products not included in the Hypothetical Monopolist to pre-merger levels.
      When TRUE, \code{defineMarkets}
      assumes that all of the products not included in the Hypothetical
      Monopolist best respond to the Hypothetical Monopolist's prices.}
     \item{startProd}{ An integer equal to the index (position) of the merging
      parties' product used to start the Hypothetical Monopolist
      Test. If missing (the default) the
      first merging parties' product listed is used.}
    \item{prodOrder}{A length-k vector of product indices whose order
      determines the sequence in which products are placed under the control
      of the Hypothetical Monopolist. The first element must always equal to
      `startProd'. Default is missing, in which case products are added to
      the Hypothetical Monopolist in order from greatest to least
      diversion from `startProd'. }
    \item{...}{Pass options to the optimizer used to solve for
      equilibrium prices. Only needed when `supplyResponse' is TRUE. }
  }

   \details{
     \code{isMarket} implements the Hypothetical Monopolist Test
     on the products indexed by `prodIndex' for a `ssnip'. The
     Hypothetical Monopolist Test determines whether a profit-maximizing
     Hypothetical Monopolist who owns the products indexed by
     `prodIndex' would increase the price of at least one of the merging
     parties' products by a
     small, significant, and non-transitory amount (i.e. a SSNIP). If it
     would, then the products indexed by `prodIndex' constitute an
     antitrust market.

     \code{smallestMarket} implements the "smallest market principle"
     for the Hypothetical Monopolist Test.  Beginning with the product
     indexed by `startProd', \code{smallestMarket} places products under
     the Hypothetical Monopolist`s control in the order dictated by
     `prodOrder' until the Hypothetical Monopolist satisfies the
     `ssnip'. \code{smallestMarket} then returns the set of products
     that satisfied the `ssnip'. When `prodOrder' is missing,
     \code{smallestMarket} adds products
      to the Hypothetical Monopolist in order from greatest to least
      diversion from `startProd'.

      By default, both \code{isMarket} and \code{smallestMarket} assume
      that the prices of all products not owned by the Hypothetical
      Monopolist are priced at their (predicted) pre-merger
      levels. Setting `supplyResponse' equal to FALSE relaxes this
      assumption, allowing their owners to set prices optimally in
      response to the Hypothetical Monopolist's prices.

      \code{calcPriceDeltaHypoMon} calculates the price changes
      relative to (predicted) pre-merger prices that a
      Hypothetical Monopolist would impose on the products indexed by
      `prodIndex',
      holding the prices of products not
      included in the Hypothetical Monopolist fixed at pre-merger
      levels. \code{calcPriceDeltaHypoMon} is called by
      \code{defineMarkets} when `supplyResponse' equals FALSE.
    }

    \value{
      \code{isMarket} returns TRUE if the products indexed by 'prodIndex'
      constitute an antitrust market under the specified 'ssnip', and
      FALSE otherwise. \code{isMarket} returns an error if 'prodIndex'
      does not contain at least one of the merging parties
      products. \code{isMarket} returns a warning if there are products
      not included in the candidate market that have higher diversion
      from  the merging parties' products included in the candidate
      market than
      those contained in 'prodIndex'.

      \code{smallestMarket} returns a vector whose first k elements
      either equal to NA if the product is not included in the smallest
      market or 1 if it is. The last element records the maximum
      percent price change for a merging parties' product placed under
      the control of the Hypothetical Monopolist.

      \code{calcPriceDeltaHypoMon} returns a vector of price changes for
      all products placed under the control of the Hypothetical
      Monopolist (i.e. all products indexed by 'prodIndex').

    }

\references{
  U.S. Department of Justice and Federal Trade Commission,
  \emph{Horizontal Merger Guidelines}. Washington DC: U.S.
  Department of Justice, 2010.
  \url{http://www.justice.gov/atr/public/guidelines/hmg-2010.html}
  (accessed July 29, 2011).
  }

