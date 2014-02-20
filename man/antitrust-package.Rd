\name{antitrust-package}
\alias{antitrust-package}
\alias{antitrust}
\docType{package}
\title{
Antitrust Library
}
\description{
  A collection of tools for antitrust practitioners, including the ability to calibrate different consumer demand systems and simulate the effects mergers under different competitive regimes.
}

\details{
\tabular{ll}{
Package: \tab antitrust\cr
Type: \tab Package\cr
Version: \tab 0.93\cr
Date: \tab 2014-02-19\cr
License: \tab Unlimited \cr
LazyLoad: \tab yes\cr
}
}

\section{Disclaimer}{
The views expressed herein are entirely those of the authors and should
not be purported to reflect those of the U.S. Department of Justice.
The \code{antitrust} package has been released into the
public domain without warranty of any kind, expressed or implied.
Address: Economic Analysis Group, Antitrust Division, U.S. Department of
Justice, 450 5th St. NW, Washington DC 20530.
E-mail: charles.taragin@usdoj.gov and michael.sandfort@usdoj.gov.
}

\section{Getting Started}{
  \enumerate{
    \item Collect data on product prices, shares, margins and diversions
    (optional).
    \item If you have data on many/all products in the market consider
    calibrating a demand system and simulating a merger with either
    a \code{\link{aids}},\code{\link{logit}}, \code{\link{ces}},
    \code{\link{linear}}, or \code{\link{loglin}} demand system.
    \item If you only have data on the merging parties' products, consider
    using \code{\link{cmcr.bertrand}} or \code{\link{cmcr.cournot}} to
    uncover the marginal cost reductions needed to offset a post-merger increase.

  }

}


\author{
Charles Taragin and Michael Sandfort

Maintainer: Charles Taragin \email{charles.taragin@usdoj.gov}
}


