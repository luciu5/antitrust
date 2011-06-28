\name{Antitrust-class}
\Rdversion{1.1}
\docType{class}
\alias{Antitrust-class}
%%\alias{cmcr,Antitrust-method}
%%\alias{diversion,Antitrust-method}
\alias{ownerToMatrix,Antitrust-method}
\alias{hhi,Antitrust-method}
\alias{show,Antitrust-method}
\alias{summary,Antitrust-method}
\alias{calcMargins,Antitrust-method}
\alias{calcPriceDelta,Antitrust-method}

\title{Class "Antitrust"}
\description{
  The "Antitrust" class is a building block used to create other classes
  in this package. As such, it is most likely to be useful for developers
  who wish to code their own merger calibration/simulation routines.


}
\section{Objects from the Class}{

Objects can be created by calls of the form \code{new("Antitrust", ...)}.

}
\section{Slots}{
  Let k denote the number of products produced by all firms.
  \describe{
    \item{\code{shares}:}{A length k vector containing quantity shares}
    \item{\code{mcDelta}:}{A length k vector where each element equals the
    proportional change in a product's marginal costs due to
    the merger.}
    \item{\code{slopes}:}{A k x (k+1) matrix of linear demand intercepts and slope coefficients}
    \item{\code{ownerPre}:}{ A k x k matrix of pre-merger ownership shares.}
    \item{\code{ownerPost}:}{ A k x k matrix of post-merger ownership shares.}
    \item{\code{labels}:}{A length k vector of labels.}

  }
}
\section{Methods}{
  For all of methods containing the 'preMerger' argument, 'preMerger' takes
  on a value of TRUE or FALSE, where TRUE invokes the method using the
  pre-merger ownership structure, while FALSE invokes the method using
  the post-merger ownership structure.
  \describe{
    \item{calcMargins}{\code{signature(object = "Antitrust", preMerger
	= TRUE)}: Compute either pre-merger or post-merger equilibrium
      margins under the assumption that firms play a differentiated product
      Bertrand Nash pricing game.}
    \item{\code{\link{calcPriceDelta}}}{\code{signature(object =
	"Antitrust")}: Calculates the proportional change in product
      prices from a merger.}
    \item{\code{\link{cmcr}}}{\code{signature(object =
	"Antitrust")}: Calculates compensated marginal cost reduction,
      the percentage decrease in the marginal costs of the merging
      parties' products needed to offset a post-merger price increase.}
     \item{hhi}{\code{signature(object = "Antitrust", preMerger
	= TRUE)}: Compute either the pre-merger or post-merger
       Herfindahl–Hirschman Index (HHI) under the assumption that firms play a differentiated product
      Bertrand Nash pricing game.}
    \item{\code{\link{diversion}}}{\code{signature(object = "Antitrust", preMerger
      = TRUE)}: Computes a k x k matrix of diversion ratios.}
    \item{ownerToMatrix}{\code{signature(object = "Antitrust", preMerger
      = TRUE)}: Converts an ownership vector (or factor) to a k x k matrix of
  1s and 0s.}
    \item{show}{\code{signature(object = "Antitrust")}: Displays
      the percentage change in prices due to the merger. }
    \item{summary}{\code{signature(object = "Antitrust")}: Summarizes
      the effect of the merger, including price and quantity changes. }
	 }
}
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
