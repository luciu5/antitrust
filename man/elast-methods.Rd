\name{elast-methods}
\docType{methods}
\alias{elast-methods}
\alias{elast}
\alias{elast,ANY-method}
\alias{elast,AIDS-method}
\alias{elast,CES-method}
\alias{elast,CESNests-method}
\alias{elast,Linear-method}
\alias{elast,LogLin-method}
\alias{elast,Logit-method}
\alias{elast,LogitNests-method}
\alias{elast,Cournot-method}

\title{Methods For Calculating Own and Cross-Price Elasticities}
\description{
 Calculate the own and cross-price elasticity between any two products in the market.
}
\usage{
\S4method{elast}{ANY}(object,preMerger=TRUE,market=FALSE)
}
\arguments{
  \item{object}{An instance of one of the classes listed above.}
  \item{preMerger}{If TRUE, calculates pre-merger price elasticities. If
    FALSE, calculates post-merger price elasticities. Default is TRUE.}
  \item{market}{If TRUE, calculates the market (aggregate) elasticity. If
    FALSE, calculates matrix of own- and cross-price elasticities. Default is FALSE.}
}

\details{

  When \sQuote{market} is FALSE, this method computes the matrix
  of own and cross-price elasticities. Element i,j of this matrix is
  the percentage change in the demand for good i from a small change in
  the price of good j. When \sQuote{market} is TRUE, this method computes the
  market (aggregate) elasticities using share-weighted prices.

  When \sQuote{preMerger} is TRUE, elasticities are
  calculated at pre-merger equilibrium prices and shares, and when \sQuote{preMerger} is FALSE, they
  are calculated at post-merger equilibrium prices and shares.
}

\value{returns a k x k matrix of own- and cross-price elasticities,
  where k is the number of products in the market}
\keyword{methods}

