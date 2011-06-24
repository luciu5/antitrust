\name{Logit-class}
\Rdversion{1.1}
\docType{class}
\alias{Logit-class}
\alias{calcPrices,Logit-method}
\alias{calcShares,Logit-method}
\alias{calcSlopes,Logit-method}
%\alias{CV,Logit-method}
%\alias{elast,Logit-method}

\title{Class "Logit"}
\description{
   The "Logit" class contains all the information needed to calibrate a Logit
  demand system and perform a merger analysis under the assumption that
  firms are playing a differentiated products Bertrand pricing game.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by using the constructor function \code{\link{logit}}.
%%  ~~ describe objects here ~~
}
\section{Slots}{
  Let k denote the number of products produced by all firms.
  \describe{
    \item{\code{prices}:}{A length k vector of product prices.}
    \item{\code{margins}:}{A length k vector of product margins, some of which may
      equal NA.}
    \item{\code{mc}:}{A length k vector containing (constant) marginal costs. }
    \item{\code{pricePre}:}{A length k vector of simulated pre-merger prices.}
    \item{\code{pricePost}:}{A length k vector of simulated post-merger prices.}
    \item{\code{priceStart}:}{A length k vector of starting values used to solve for
    equilibrium price.}
    \item{\code{normIndex}:}{An integer specifying the product index against which
     the mean values of all other products are normalized.}
    \item{\code{shareInside}:}{The share of customers that purchase any of the
      products included in the 'prices' vector.}
     \item{\code{slopes}:}{A list containing the coefficient on price ('alpha') and the vector of mean
      valuations ('meanval')}

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
    \item{calcPrices}{\code{signature(object = "Logit", preMerger
	= TRUE)}: Compute either pre-merger or post-merger equilibrium
      prices under the assumptions that consumer demand is Logit and firms play a differentiated product
      Bertrand Nash pricing game.}
    \item{calcShares}{\code{signature(object = "Logit", preMerger
	= TRUE)}: Compute either pre-merger or post-merger equilibrium
      quantity shares under the assumptions that consumer demand is Logit and firms play a differentiated product
      Bertrand Nash pricing game.}
    \item{calcSlopes}{\code{signature(object = "Logit")}: Uncover Logit demand
      parameters. Assumes that firms are currently at equilibrium  in a
      differentiated product Bertrand Nash pricing game.}
    \item{\code{\link{CV}}}{\code{signature(object =
	"Logit")}: Calculate the amount of money a representative consumer would need to
      be paid to be just as well off as they were before the merger.}
    \item{\code{\link{elast}}}{\code{signature(object = "Logit", preMerger
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
showClass("Logit")
}
\keyword{classes}
