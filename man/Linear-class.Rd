\name{Linear-class}
\Rdversion{1.1}
\docType{class}
\alias{Linear-class}
%\alias{calcPrices,Linear-method}
%\alias{calcQuantities,Linear-method}
%\alias{calcShares,Linear-method}
%\alias{calcSlopes,Linear-method}

\title{Class \dQuote{Linear}}
\description{
   The \dQuote{Linear} class contains all the information needed to
  calibrate a Linear
  demand system and perform a merger simulation analysis under the assumption that
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
     \item{\code{intercepts}:}{A length k vector of demand intercepts.}
     \item{\code{prices}:}{A length k vector product prices.}
     \item{\code{quantities}:}{A length k vector of product quantities.}
     \item{\code{margins}:}{A length k vector of product margins. All margins must
       be between 0 and 1.}
     \item{\code{diversion}:}{A k x k matrix of diversion ratios with diagonal
       elements equal to 1.}
     \item{\code{priceStart}:}{A length k vector of prices used as the initial guess
       in the nonlinear equation solver.}
     \item{\code{symmetry}:}{If TRUE, requires the matrix of demand slope coefficients to be
       consistent with utility maximization theory.}

   }
 }
 \section{Extends}{
   Class \code{\linkS4class{Bertrand}}, directly.
   Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 2.
 }
 \section{Methods}{
   For all of methods containing the \sQuote{preMerger} argument, \sQuote{preMerger} takes
   on a value of TRUE or FALSE, where TRUE invokes the method using the
   pre-merger ownership structure, while FALSE invokes the method using
   the post-merger ownership structure.
   \describe{
     \item{\code{calcPrices}}{\code{signature(object, preMerger
	 = TRUE,...)}}{ Compute either pre-merger or post-merger equilibrium
       prices under the assumptions that consumer demand is Logit and firms play a differentiated product
       Bertrand Nash pricing game. \sQuote{...} may be used to change the
       default values of \code{\link{constrOptim}}, the non-linear
       equation solver used to enforce non-negative equilibrium quantities. }
     \item{\code{\link{calcPriceDeltaHypoMon}}}{\code{signature(object,prodIndex,...)}}{ Calculates the price changes that a
       Hypothetical Monopolist would impose on its products
       relative to pre-merger prices. }
     \item{\code{calcQuantities}}{\code{signature(object, preMerger
	 = TRUE)}}{ Compute either pre-merger or post-merger equilibrium
       quantities under the assumptions that consumer demand is Linear and firms play a differentiated product
       Bertrand Nash pricing game.}
     \item{\code{calcShares}}{\code{signature(object, preMerger
	 = TRUE, revenue = FALSE)}}{ Compute either pre-merger or post-merger equilibrium
       quantity shares under the assumptions that consumer demand is Linear and firms play a differentiated product
       Bertrand Nash pricing game.}
     \item{\code{calcSlopes}}{\code{signature(object)}}{ Uncover
       slopes and intercept from a Linear demand system. Assumes that firms are currently at equilibrium  in a
       differentiated product Bertrand Nash pricing game.}
     \item{\code{\link{CV}}}{\code{signature(object =
	 "Linear")}}{ Calculate the amount of money a representative consumer would need to
       be paid to be just as well off as they were before the merger.}
     \item{\code{\link{elast}}}{\code{signature(object, preMerger
	 = TRUE)}}{ Computes a k x k matrix of own and cross-price elasticities.}
   }
 }

\author{Charles Taragin \email{charles.taragin@usdoj.gov}}

\examples{
showClass("Linear")           # get a detailed description of the class
showMethods(classes="Linear") # show all methods defined for the class
}
\keyword{classes}
