% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AntitrustClasses.R
\name{Antitrust-Class}
\alias{Antitrust-Class}
\alias{Antitrust-class}
\alias{matrixOrList-class}
\alias{matrixOrVector-class}
\alias{characterOrList-class}
\title{\dQuote{Antitrust} Classes}
\description{
The \dQuote{Antitrust} class is a building block used to create other classes
in this package. As such, it is most likely to be useful for developers
who wish to code their own calibration/simulation routines.

Let k denote the number of products produced by all firms below.
}
\section{Slots}{

\describe{
\item{\code{pricePre}}{A length k vector of simulated pre-merger prices.}

\item{\code{pricePost}}{A length k vector of simulated post-merger prices.}

\item{\code{output}}{A length 1 logical vector equal to TRUE if merger simulates the effect in an output market,
and FALSE if the merger simulates the effect in an input market. Default is TRUE.}

\item{\code{ownerPre}}{A k x k matrix of pre-merger ownership shares.}

\item{\code{ownerPost}}{A k x k matrix of post-merger ownership shares.}

\item{\code{labels}}{A length k vector of labels.}

\item{\code{control.slopes}}{A list of \code{\link{optim}} control parameters passed to the calibration routine optimizer (typically the \code{calcSlopes} method).}

\item{\code{control.equ}}{A list of \code{\link[BB]{BBsolve}} control parameters passed to the non-linear equation solver (typically the \code{calcPrices} method).}
}}

\section{Objects from the Class}{

Objects can be created by calls of the form \code{new("Antitrust", ...)}.
}

\section{The \dQuote{matrixOrList}, \dQuote{matrixOrVector} and \dQuote{characterOrList} Classes}{

The \dQuote{matrixOrList},\dQuote{matrixOrVector} and \dQuote{characterOrList} classes are
virtual classes used for validity checking in the \sQuote{ownerPre} and
\sQuote{ownerPost} slots of \dQuote{Antitrust} and the \sQuote{slopes} slot in
\dQuote{Bertrand}.
}

\examples{
showClass("Antitrust")           # get a detailed description of the class
}
\author{
Charles Taragin \email{ctaragin+antitrustr@gmail.com}
}
\keyword{classes}
