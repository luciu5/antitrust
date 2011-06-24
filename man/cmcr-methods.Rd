\name{cmcr-methods}
\docType{methods}
\alias{cmcr-methods}
\alias{cmcr}
\alias{cmcr,Antitrust-method}
\alias{cmcr,PCAIDS-method}
\title{Methods For Calculating Compensating Marginal Cost Reductions (Bertrand)}

\description{
  Calculate the marginal cost reductions necessary to restore
  premerger prices in a merger involving firms playing a differentiated
  products Bertrand pricing game.

}
\section{Methods}{

\describe{

\item{\code{signature(object = c("Antitrust","PCAIDS"))}}{
This method uses the results from the merger simulation and calibration
methods associates with a particular class to compute the compensating
marginal cost reduction (CMCR).
}
}}


\seealso{\code{\link{cmcr.bertrand}} is a function that calculates CMCR
  without the need to first calibrate a demand system and simulate a merger.}

\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
