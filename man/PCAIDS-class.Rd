\name{PCAIDS-class}
\Rdversion{1.1}
\docType{class}
\alias{PCAIDS-class}
\alias{calcMargins,PCAIDS-method}
\alias{calcPriceDelta,PCAIDS-method}
\alias{calcShares,PCAIDS-method}
\alias{calcSlopes,PCAIDS-method}
%%\alias{cmcr,PCAIDS-method}
%%\alias{CV,PCAIDS-method}
%%\alias{elast,PCAIDS-method}
\alias{show,PCAIDS-method}
\alias{summary,PCAIDS-method}

\title{Class "PCAIDS"}
\description{
  The "PCAIDS" class contains all the information needed to
  calibrate a PCAIDS
  demand system and perform a merger analysis under the assumption that
  firms are playing a differentiated products Bertrand pricing game.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by using the constructor function \code{\link{pcaids}}.
%%  ~~ describe objects here ~~
}
\section{Slots}{
   Let k denote the number of products produced by all firms.
  \describe{
     \item{\code{knownElast}:}{A negative number equal to the pre-merger own-price
       elasticity for any of the k products.}
     \item{\code{knownElastIndex}:}{An integer equal to the position of the
    'knownElast' product in the 'shares' vector.}
    \item{\code{mktElast}:}{A negative number equal to the industry pre-merger
      price elasticity.}
    \item{\code{priceDeltaStart}:}{A length k vector who elements equal to an
      initial guess of the proportional change in price caused by the
      merger.}
    \item{\code{priceDelta}:}{A length k vector containing the
      simulated price effects from the merger.}
    \item{\code{diversion}:}{A k x k matrix  containing product diversions.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{Antitrust}"}, directly.
}
\section{Methods}{
  For all of methods containing the 'preMerger' argument, 'preMerger' takes
  on a value of TRUE or FALSE, where TRUE invokes the method using the
  pre-merger ownership structure, while FALSE invokes the method using
  the post-merger ownership structure.
  \describe{
    \item{calcMargins}{\code{signature(object = "PCAIDS",
	preMerger=TRUE)}: Calculates pre-merger or post-merger
      equilibrium margins. }
    \item{calcPriceDelta}{\code{signature(object = "PCAIDS")}:
      Computes the proportional change in each products' price from the
      merger under the assumptions that consumer demand is PCAIDS and firms play a differentiated product
      Bertrand Nash pricing game.}
    \item{calcShares}{\code{signature(object = "PCAIDS", preMerger
	= TRUE)}: Computes either pre-merger or post-merger equilibrium
      quantity shares under the assumptions that consumer demand is PCAIDS and firms play a differentiated product
      Bertrand Nash pricing game.}
    \item{calcSlopes}{\code{signature(object = "PCAIDS")}: Uncover PCAIDS demand
      parameters. Assumes that firms are currently at equilibrium  in a
      differentiated product Bertrand Nash pricing game.}
    \item{\code{\link{cmcr}}}{\code{signature(object = "PCAIDS")}: Calculates compensated marginal cost reduction,
      the percentage decrease in the marginal costs of the merging
      parties' products needed to offset a post-merger price increase. }
    \item{\code{\link{CV}}}{\code{signature(object =
	"PCAIDS", prices)}: Calculate the amount of money a representative consumer would need to
      be paid to be just as well off as they were before the
      merger. Requires a length-k vector of pre-merger prices.}
    \item{\code{\link{elast}}}{\code{signature(object = "PCAIDS", preMerger
	= TRUE)}: Computes a k x k matrix of own and cross-price
      elasticities.}
    \item{show}{\code{signature(object = "PCAIDS")}: Displays
      the percentage change in prices due to the merger. }
    \item{summary}{\code{signature(object = "PCAIDS")}: Summarizes
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
showClass("PCAIDS")
}
\keyword{classes}
