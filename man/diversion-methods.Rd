\name{diversion-methods}
\docType{methods}
\alias{diversion-methods}
\alias{diversion}
\alias{diversion,Antitrust-method}
\alias{diversion,CES-method}
\title{Methods For Calculating Diversion}
\description{
 Calculate the diversion matrix between any two products in the market.
}
\section{Methods}{
\describe{

\item{\code{signature(object = c("Antitrust","CES"), preMerger=TRUE)}}{
  %%  ~~describe this method her
  This method uses the results from the merger calibration and
  simulation to compute the diversion matrix between any two products
  in the market. Element i,j of this matrix is the diversion from
  product i to product j, or the
  proportion of product i's sales that leave (go to) i for (from) j due
  to a increase (decrease) in i's price. Mathematically, diversion is
  \eqn{\frac{-\epsilon_{ji}share_j}{\epsilon_{ii}share_i}}{(-e_ji/e_ii)(share_j/share_i)},
  where \eqn{\epsilon_{ji}}{e_ji} is the cross-price elasticity from j
  to i. When 'preMerger' is TRUE, diversions are
  calculated at pre-merger equilbrium prices, and when 'preMerger' is FALSE, they
  are calculated at post-merger equilibrium prices.
}
}
}

\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
