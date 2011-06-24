\name{LogLog-class}
\Rdversion{1.1}
\docType{class}
\alias{LogLog-class}
\alias{calcPrices,LogLog-method}
\alias{calcQuantities,LogLog-method}
\alias{calcSlopes,LogLog-method}
%%\alias{elast,LogLog-method}

\title{Class "LogLog"}
\description{
  The "LogLog" class contains all the information needed to
  calibrate a Log-Linear
  demand system and perform a merger analysis under the assumption that
  firms are playing a differentiated Bertrand products pricing game.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by using the constructor function \code{\link{loglog}}.
%%  ~~ describe objects here ~~
}
\section{Slots}{
   Let k denote the number of products produced by all firms.
  \describe{
    \item{\code{priceStart}:}{A length k vector of prices used as the initial guess
  in the nonlinear equation solver.}

  }
}
\section{Extends}{
Class \code{"\linkS4class{Linear}"}, directly.
Class \code{"\linkS4class{Antitrust}"}, by class "Linear", distance 2.
}
\section{Methods}{
  For all of methods containing the 'preMerger' argument, 'preMerger' takes
  on a value of TRUE or FALSE, where TRUE invokes the method using the
  pre-merger ownership structure, while FALSE invokes the method using
  the post-merger ownership structure.
  \describe{
    \item{calcPrices}{\code{signature(object = "LogLog", preMerger
	= TRUE)}: Compute either pre-merger or post-merger equilibrium
      prices under the assumptions that consumer demand is Log-Linear and firms play a differentiated product
      Bertrand Nash pricing game.}
    \item{calcQuantities}{\code{signature(object = "LogLog", preMerger
	= TRUE)}: Compute either pre-merger or post-merger equilibrium
      quantities under the assumptions that consumer demand is Log-Linear and firms play a differentiated product
      Bertrand Nash pricing game.}
    \item{calcSlopes}{\code{signature(object = "LogLog")}: Uncover
      slopes and intercept from a Log-Linear demand system. Assumes that firms are currently at equilibrium  in a
      differentiated product Bertrand Nash pricing game.}
    \item{\code{\link{elast}}}{\code{signature(object = "LogLog", preMerger
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
showClass("LogLog")
}
\keyword{classes}
