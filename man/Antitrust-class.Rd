\name{Antitrust-class}
\Rdversion{1.1}
\docType{class}
\alias{Antitrust-class}
\alias{matrixOrList-class}
\alias{matrixOrVector-class}
\alias{characterOrList-class}
%\alias{cmcr,Antitrust-method}
%\alias{diversion,Antitrust-method}
%\alias{hhi,Antitrust-method}
\alias{show,Antitrust-method}
%\alias{summary,Antitrust-method}


\title{Class \dQuote{Antitrust}}
\description{
  The \dQuote{Antitrust} class is a building block used to create other classes
  in this package. As such, it is most likely to be useful for developers
  who wish to code their own calibration/simulation routines.


}
\section{Objects from the Class}{

Objects can be created by calls of the form \code{new("Antitrust", ...)}.

}
\section{Slots}{
  Let k denote the number of products produced by all firms.
  \describe{
    \item{\code{pricePre}:}{A length k vector of simulated pre-merger prices.}
    \item{\code{pricePost}:}{A length k vector of simulated post-merger prices.}
    \item{\code{ownerPre}:}{ A k x k matrix of pre-merger ownership shares.}
    \item{\code{ownerPost}:}{ A k x k matrix of post-merger ownership shares.}
    \item{\code{labels}:}{A length k vector of labels.}
    \item{\code{control.slopes}:}{A list of  \code{\link{optim}}  control parameters passed to the calibration routine optimizer (typically the \code{calcSlopes} method).}
    \item{\code{control.equ}:}{A list of  \code{\link[BB]{BBsolve}} control parameters passed to the non-linear equation solver (typically the \code{calcPrices} method).}
  }
}
\section{Methods}{
  Many of the methods described below contain a \sQuote{preMerger} argument. The \sQuote{preMerger} takes
  on a value of TRUE or FALSE, where TRUE invokes the method using the
  pre-merger values, while FALSE invokes the method using
  the post-merger ownership structure.
  \describe{
    \item{\code{calcPriceDelta}}{\code{signature(object)}}{ Calculates the proportional change in product
      prices from a merger.}
    \item{\code{ownerToMatrix}}{\code{signature(object, preMerger
	= TRUE)}}{ Converts an ownership vector (or factor) to a k x k matrix of
      1s and 0s.}
     \item{\code{ownerToVec}}{\code{signature(object, preMerger
	= TRUE)}}{ Converts a k x k  ownership matrix to a length-k
      vector whose values identify an owner.}
    \item{\code{show}}{\code{signature(object)}}{ Displays
      the percentage change in prices due to the merger. }
     }}

\section{The \dQuote{matrixOrList}, \dQuote{matrixOrVector} and \dQuote{characterOrList} Classes}{
The \dQuote{matrixOrList},\dQuote{matrixOrVector} and \dQuote{characterOrList} classes are
virtual classes used for validity checking in the \sQuote{ownerPre} and
\sQuote{ownerPost} slots of \dQuote{Antitrust} and the \sQuote{slopes} slot in
\dQuote{Bertrand}.
}

\author{Charles Taragin \email{charles.taragin@usdoj.gov}}


\examples{
showClass("Antitrust")           # get a detailed description of the class
showMethods(classes="Antitrust") # show all methods defined for the class
}
\keyword{classes}
