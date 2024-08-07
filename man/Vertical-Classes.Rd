% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VerticalClasses.R
\docType{class}
\name{Vertical-Classes}
\alias{Vertical-Classes}
\alias{VertBargBertLogit}
\alias{VertBarg2ndLogit}
\alias{VertBargBertLogitNests}
\alias{VertBarg2ndLogitNests}
\alias{VertBargBertLogit-class}
\alias{VertBargBertLogitNests-class}
\alias{VertBarg2ndLogit-class}
\alias{VertBarg2ndLogitNests-class}
\title{\dQuote{Vertical} Classes}
\description{
The \dQuote{Vertical} classes are building blocks used to create other classes
in this package. As such, it is most likely to be useful for developers
who wish to code their own calibration/simulation routines.

The \dQuote{VertBargBertLogit} class has the information for a Vertical Supply Chain with Logit demand and a downstream Nash-Bertrand Pricing game.

The \dQuote{VertBarg2ndLogit} class has the information fora Vertical Supply Chain with Logit demand and a downstream 2nd Score Auction.
}
\section{Slots}{

\describe{
\item{\code{up}}{an instance of \dQuote{Bargaining} class.}

\item{\code{down}}{For  \dQuote{VertBargBertLogit}, an instance of  \dQuote{Logit} class.For  \dQuote{VertBarg2ndLogit}, an instance of  \dQuote{Auction2ndLogit} class.}

\item{\code{constrain}}{A length 1 character vector equal to "global", "pair", "wholesaler", or "retailer.}

\item{\code{chain_level}}{A length 1 character vector equal to "full", "wholesaler", or "retailer".}
}}

\author{
Charles Taragin \email{ctaragin+antitrustr@gmail.com}
}
\keyword{classes}
