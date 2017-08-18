\name{AIDS-class}
\Rdversion{1.1}
\docType{class}
\alias{AIDS-class}

%\alias{calcMargins,AIDS-method}
%\alias{calcPrices,AIDS-method}
%\alias{calcPriceDelta,AIDS-method}
%\alias{calcShares,AIDS-method}
%\alias{calcSlopes,AIDS-method}
%\alias{CV,AIDS-method}
%\alias{elast,AIDS-method}
\alias{show,AIDS-method}
%\alias{summary,AIDS-method}

\title{Class \dQuote{AIDS}}
\description{
  The \dQuote{AIDS} class contains all the information needed to
  calibrate a AIDS
  demand system and perform a merger analysis under the assumption that
  firms are playing a differentiated products Bertrand pricing game.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by using the constructor function \code{\link{aids}}.
%%  ~~ describe objects here ~~
}
\section{Slots}{
   Let k denote the number of products produced by all firms.
  \describe{
    \item{\code{mktElast}:}{A negative number equal to the industry pre-merger
      price elasticity.}
    \item{\code{priceStart}:}{A length k vector who elements equal to an
      initial guess of the proportional change in prices caused by the
      merger.}
      \item{\code{parmStart}:}{A length 2 vector who elements equal to an
      initial of a single diagonal element of the matrix of slope coefficients, as well as the market elasticity. }
    \item{\code{priceDelta}:}{A length k vector containing the
      simulated price effects from the merger.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{Linear}}, directly.
  Class \code{\linkS4class{Bertrand}}, by class \dQuote{Linear}, distance 2.
}
\section{Methods}{
  For all of methods containing the \sQuote{preMerger} argument, \sQuote{preMerger} takes
  on a value of TRUE or FALSE, where TRUE invokes the method using the
  pre-merger ownership structure, while FALSE invokes the method using
  the post-merger ownership structure.
  \describe{
    \item{\code{calcMargins}}{\code{signature(object ,
	preMerger=TRUE)}}{ Calculates pre-merger or post-merger
      equilibrium margins. }
    \item{\code{calcPriceDelta}}{\code{signature(object,isMax=FALSE,...)}}{
      Computes the proportional change in each products' price from the
      merger under the assumptions that consumer demand is AIDS and firms play a differentiated product
      Bertrand Nash pricing game.When isMax equals TRUE, a check is run
      to determine if the calculated equilibrium price vector locally maximizes
      profits. \sQuote{...} may be used to change the
      default values of \code{\link[BB]{BBsolve}}, the non-linear equation solver.}
    \item{\code{calcPrices}}{\code{signature(object, preMerger
	= TRUE)}}{ Compute either pre-merger or post-merger equilibrium
      prices under the assumptions that consumer demand is AIDS and firms play a differentiated product
      Bertrand Nash pricing game. return a vector of length-k vector of
      NAs if user did not supply prices. }
    \item{\code{\link{calcPriceDeltaHypoMon}}}{\code{signature(object,prodIndex,...)}}{ Calculates the price changes that a
   Hypothetical Monopolist would impose on its products
   relative to pre-merger prices.}
    \item{calcShares}{\code{signature(object, preMerger
	= TRUE)}}{ Computes either pre-merger or post-merger equilibrium
      quantity shares under the assumptions that consumer demand is AIDS and firms play a differentiated product
      Bertrand Nash pricing game.}
    \item{\code{calcSlopes}}{\code{signature(object)}}{ Uncover AIDS demand
      parameters. Assumes that firms are currently at equilibrium  in a
      differentiated product Bertrand Nash pricing game.}
    \item{\code{\link{cmcr}}}{\code{signature(object)}}{ Calculates compensated marginal cost reduction,
      the percentage decrease in the marginal costs of the merging
      parties' products needed to offset a post-merger price increase. }
    \item{\code{\link{CV}}}{\code{signature(object)}}{ Calculate the amount of money a representative consumer would need to
      be paid to be just as well off as they were before the
      merger. Requires a length-k vector of pre-merger prices.}
    \item{\code{\link{diversion}}}{\code{signature(object, preMerger= TRUE)}}{ Computes a k x k matrix of diversion ratios.}
    \item{\code{\link{elast}}}{\code{signature(object , preMerger
	= TRUE)}}{ Computes a k x k matrix of own and cross-price
      elasticities.}
    \item{\code{show}}{\code{signature(object)}}{ Displays
      the percentage change in prices due to the merger. }
    \item{summary}{\code{signature(object,revenue=TRUE,parameters=FALSE,digits=2,..)}}{ Summarizes
      the effect of the merger, including price and revenue
      changes. Setting \sQuote{revenue} equal to FALSE reports quantity rather
      than revenue changes. Setting \sQuote{parameters} equal to TRUE reports
      all demand parameters. \sQuote{digits} controls the number of significant
      digits reported in output.}
  }
}


\author{Charles Taragin \email{charles.taragin@usdoj.gov}}

\examples{
showClass("AIDS")           # get a detailed description of the class
showMethods(classes="AIDS") # show all methods defined for the class
}
\keyword{classes}
