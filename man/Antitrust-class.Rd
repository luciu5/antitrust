\name{Antitrust-class}
\Rdversion{1.1}
\docType{class}
\alias{Antitrust-class}
%\alias{cmcr,Antitrust-method}
%\alias{diversion,Antitrust-method}
%\alias{hhi,Antitrust-method}
\alias{show,Antitrust-method}
%\alias{summary,Antitrust-method}


\title{Class "Antitrust"}
\description{
  The "Antitrust" class is a building block used to create other classes
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

  }
}
\section{Methods}{
  Many of the methods described below contain a `preMerger' argument. The `preMerger' takes
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
\references{
%%  ~~put references to the literature/web site here~~
}
\author{Charles Taragin \email{charles.taragin@usdoj.gov}}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("Antitrust")
}
\keyword{classes}
