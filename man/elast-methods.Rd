\name{elast-methods}
\docType{methods}
\alias{elast-methods}
\alias{elast}
\alias{elast,CES-method}
\alias{elast,CESNests-method}
\alias{elast,Linear-method}
\alias{elast,Logit-method}
\alias{elast,LogitNests-method}
\alias{elast,LogLog-method}
\alias{elast,PCAIDS-method}
\title{Methods For Calculating Own and Cross-Price Elasticities}
\description{
 Calculate the own and cross-price elasticity between any two products in the market.
}
\section{Methods}{
\describe{

\item{\code{signature(object = c("CES","CESNests",
    "Linear","Logit","LogitNests","LogLog","PCAIDS"), preMerger=TRUE)}}{
This method uses the results from the merger to compute the matrix
of own and cross-price elasticities. Element i,j of this matrix is
the percentage change in the demand for good i from a small change in
    the price of good j. When 'preMerger' is TRUE, elasticities are
    calculated at pre-merger equilbrium prices, and when 'preMerger' is FALSE, they
    are calculated at post-merger equilibrium prices.
}
}}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
