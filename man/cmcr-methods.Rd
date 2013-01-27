\name{cmcr-methods}
\docType{methods}
\alias{cmcr-methods}
\alias{upp-methods}
\alias{cmcr,ANY-method}
\alias{cmcr,Bertrand-method}
\alias{cmcr,AIDS-method}
\alias{upp,ANY-method}
\alias{upp,Bertrand-method}
\alias{upp,AIDS-method}

\title{Methods For Calculating Compensating Marginal Cost Reductions and
Upwards Pricing Pressure Index (Bertrand)}

\description{
  Calculate the marginal cost reductions necessary to restore
  premerger prices in a merger, or the Upwards Pricing Pressure Index
  for the products of merging firms playing a differentiated
  products Bertrand pricing game.
}

\usage{
\S4method{cmcr}{ANY}(object)
\S4method{upp}{ANY}(object)
}

\arguments{
  \item{object}{An instance of one of the classes listed above.}
}

\details{
  \code{cmcr} uses the results from the merger simulation and calibration
    methods associates with a particular class to compute the compensating
    marginal cost reduction (CMCR) for each of the merging parties' products.

    Like \code{cmcr}, \code{upp} uses the results from the merger simulation and
    calibration to compute the upwards pricing pressure of the merger on
    each merging parties' products.

}

\value{
  \code{cmcr} returns a vector of length k equal to CMCR for the
  merging parties' products and 0 for all other products.

  \code{upp} returns a vector of length k equal to the  net UPP for the
  merging parties' products and 0 for all other products.
}

\seealso{\code{\link{cmcr.bertrand}} is a function that calculates CMCR
  without the need to first calibrate a demand system and simulate a
  merger. Likewise,\code{\link{upp.bertrand}} calculates net UPP
  without the need to first calibrate a demand
  system and simulate a merger.}

\keyword{methods}

