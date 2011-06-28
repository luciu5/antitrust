\name{Linear-class}
\Rdversion{1.1}
\docType{class}
\alias{Linear-class}
\alias{calcPrices,Linear-method}
\alias{calcQuantities,Linear-method}
\alias{calcShares,Linear-method}
\alias{calcSlopes,Linear-method}
%%\alias{CV,Linear-method}
%%\alias{elast,Linear-method}

\title{Class "Linear"}
\description{
   The "Linear" class contains all the information needed to
  calibrate a Linear
  demand system and perform a merger analysis under the assumption that
  firms are playing a differentiated Bertrand products pricing game.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by using the constructor function \code{\link{linear}}.
%%  ~~ describe objects here ~~
}
\section{Slots}{
   Let k denote the number of products produced by all firms.
  \describe{
    \item{prices}{A length k vector product prices.}
    \item{quantities}{A length k vector of product quantities.}
    \item{margins}{A length k vector of product margins. All margins must
      be between 0 and 1.}
    \item{\code{mc}:}{A length k vector containing (constant) marginal costs. }
    \item{\code{pricePre}:}{A length k vector of simulated pre-merger prices.}
    \item{\code{pricePost}:}{A length k vector of simulated post-merger prices.}
    \item{\code{diversion}:}{A k x k matrix of diversion ratios with diagonal
      elements equal to 1.}
     \item{\code{symmetry}:}{If TRUE, requires the matrix of demand slope coefficients to be
    consistent with utility maximization theory.}

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
    \item{calcPrices}{\code{signature(object = "Linear", preMerger
	= TRUE)}: Compute either pre-merger or post-merger equilibrium
      prices under the assumptions that consumer demand is Linear and firms play a differentiated product
      Bertrand Nash pricing game.}
    \item{calcQuantities}{\code{signature(object = "Linear", preMerger
	= TRUE)}: Compute either pre-merger or post-merger equilibrium
      quantities under the assumptions that consumer demand is Linear and firms play a differentiated product
      Bertrand Nash pricing game.}
    \item{calcShares}{\code{signature(object = "Linear", preMerger
	= TRUE)}: Compute either pre-merger or post-merger equilibrium
      quantity shares under the assumptions that consumer demand is Linear and firms play a differentiated product
      Bertrand Nash pricing game.}
    \item{calcSlopes}{\code{signature(object = "Linear")}: Uncover
      slopes and intercept from a Linear demand system. Assumes that firms are currently at equilibrium  in a
      differentiated product Bertrand Nash pricing game.}
    \item{\code{\link{CV}}}{\code{signature(object =
	"Linear")}: Calculate the amount of money a representative consumer would need to
      be paid to be just as well off as they were before the merger.}
    \item{\code{\link{elast}}}{\code{signature(object = "Linear", preMerger
	= TRUE)}: Computes a k x k matrix of own and cross-price elasticities.}
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
showClass("Linear")
}
\keyword{classes}
