\name{diversion-methods}
\docType{methods}
\alias{diversion-methods}
\alias{diversion}
\alias{diversion,ANY-method}
\alias{diversion,AIDS-method}
\alias{diversion,Bertrand-method}
\title{Methods For Calculating Diversion}
\description{
 Calculate the diversion matrix between any two products in the market.
}
\usage{
\S4method{diversion}{ANY}(object,preMerger=TRUE,revenue=FALSE)
}
\arguments{
  \item{object}{An instance of one of the classes listed above.}
  \item{preMerger}{If TRUE, calculates pre-merger price elasticities. If
    FALSE, calculates post-merger price elasticities. Default is TRUE.}
  \item{revenue}{If TRUE, calculates revenue diversion. If
    FALSE, calculates quantity diversion. Default is TRUE for `Bertrand'
    and FALSE for `AIDS'.}
}


\section{Methods}{
\describe{
  \item{diversion}{\code{signature(object=Bertrand,preMerger=TRUE,revenue=FALSE)}}{
    When `revenue' is FALSE (the default), this method uses the results from the merger calibration and
    simulation to compute the \emph{quantity} diversion matrix between any two products
    in the market. Element i,j of this matrix is the quantity diversion from
    product i to product j, or the
    proportion of product i's sales that leave (go to) i for (from) j due
    to a increase (decrease) in i's price. Mathematically, quantity diversion is
    \eqn{\frac{-\epsilon_{ji}share_j}{\epsilon_{ii}share_i}},
    where \eqn{\epsilon_{ij}} is the cross-price elasticity from i
    to j.

    When `revenue' is TRUE, this method computes the revenue diversion
    matrix between any two products in the market. Element i,j of this matrix is the revenue diversion from
    product i to product j, or the
    proportion of product i's revenues that leave (go to) i for (from) j due
    to a increase (decrease) in i's price. Mathematically, revenue
    diversion is
    \eqn{-\frac{\epsilon_{ji}(\epsilon_{jj}-1)r_j}{\epsilon_{jj}(\epsilon_{ii}-1)r_j}}
    where \eqn{r_i} is the revenue share of product i.

    When `preMerger' is TRUE, diversions are
    calculated at pre-merger equilibrium prices, and when `preMerger' is FALSE, they
    are calculated at post-merger equilibrium prices.
  }
  \item{diversion}{\code{signature(object=AIDS,preMerger=TRUE,revenue=TRUE)}}{
     When `revenue' is TRUE (the default), this method computes the \emph{revenue} diversion matrix between any two
    products in the market. For AIDS, the revenue diversion from i to j is
    \eqn{\frac{\beta_{ji}}{\beta_ij}}, where \eqn{\beta_{ij}} is the
    percentage change in product i's revenue due to a change in j's
    price.

    When `revenue' is FALSE, this \code{callNextMethod} is invoked. Will
    yield a matrix of NAs if the user did not supply prices.

    When `preMerger' is TRUE, diversions are
    calculated at pre-merger equilibrium prices, and when `preMerger' is FALSE, they
    are calculated at post-merger equilibrium prices.
  }
}
}

\value{returns a k x k matrix of diversion ratios, where the i,jth
  element is the diversion from i to j.}
\keyword{methods}

