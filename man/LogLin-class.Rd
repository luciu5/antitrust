\name{LogLin-class}
\Rdversion{1.1}
\docType{class}
\alias{LogLin-class}
%\alias{calcPrices,LogLin-method}
%\alias{calcQuantities,LogLin-method}
%\alias{calcSlopes,LogLin-method}


\title{Class \dQuote{LogLin}}
\description{
  The \dQuote{LogLin} class contains all the information needed to
  calibrate a Log-Linear
  demand system and perform a merger simulation analysis under the assumption that
  firms are playing a differentiated Bertrand products pricing game.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by using the constructor function \code{\link{loglin}}.
%%  ~~ describe objects here ~~
}
\section{Slots}{
  \describe{
    \item{\code{symmetry}:}{If TRUE, requires the matrix of demand slope coefficients to be
       consistent with utility maximization theory Default is FALSE}

  }
}
\section{Extends}{
Class \code{\linkS4class{Linear}}, directly.
Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Linear}}, distance 2.
Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 3.
}
\section{Methods}{
  For all of methods containing the \sQuote{preMerger}  argument, \sQuote{preMerger} takes
  on a value of TRUE or FALSE, where TRUE invokes the method using the
  pre-merger ownership structure, while FALSE invokes the method using
  the post-merger ownership structure.
  \describe{
    \item{\code{calcPrices}}{\code{signature(object, preMerger
	= TRUE)}}{ Compute either pre-merger or post-merger equilibrium
      prices under the assumptions that consumer demand is Log-Linear and firms play a differentiated product
      Bertrand Nash pricing game.}
    \item{\code{\link{calcPriceDeltaHypoMon}}}{\code{signature(object,prodIndex,...)}}{ Calculates the price changes that a
	Hypothetical Monopolist would impose on its products
	relative to pre-merger prices. }
      \item{\code{calcQuantities}}{\code{signature(object, preMerger
	= TRUE)}}{ Compute either pre-merger or post-merger equilibrium
      quantities under the assumptions that consumer demand is Log-Linear and firms play a differentiated product
      Bertrand Nash pricing game.}
    \item{\code{calcSlopes}}{\code{signature(object)}}{ Uncover
      slopes and intercept from a Log-Linear demand system. Assumes that firms are currently at equilibrium  in a
      differentiated product Bertrand Nash pricing game.}
    \item{\code{\link{elast}}}{\code{signature(object, preMerger
	= TRUE)}}{ Computes a k x k matrix of own and cross-price elasticities.}
	 }
}

\author{Charles Taragin \email{charles.taragin@usdoj.gov}}

\examples{
showClass("LogLin")           # get a detailed description of the class
showMethods(classes="LogLin") # show all methods defined for the class
}
\keyword{classes}
