\name{antitrust-package}
\alias{antitrust-package}
\alias{antitrust}
\docType{package}
\title{
Antitrust Library
}
\description{
  A collection of tools for antitrust practictioners, including the ability to calibrate different consumer demand systems and simulate the effects mergers under different competitive regimes.
}

\details{
\tabular{ll}{
Package: \tab antitrust\cr
Type: \tab Package\cr
Version: \tab 0.7\cr
Date: \tab 2011-01-05\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
}


\section{Getting Started}{
  \enumerate{
    \item Collect data on product prices, shares, margins and diversions
    (optional).
    \item If you have data on many/all products in the market consider
    calibrating a demand system and simulating a merger with either
    a \code{\link{pcaids}},\code{\link{logit}}, \code{\link{ces}}, or
    \code{\link{linear}} demand system.
    \item If you only have data on the merging parties' products, consider
    using \code{\link{cmcr.bertrand}} or \code{\link{cmcr.cournot}} to
    uncover the marginal cost reductions needed to offset a post-merger increase.

  }

}


\author{
Michael Sandfort and Charles Taragin

Maintainer: Charles Taragin <charles.taragin@usdoj.gov>

}
\references{
}

\seealso{
}

