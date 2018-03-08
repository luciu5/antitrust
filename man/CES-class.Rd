\name{CES-class}
\Rdversion{1.1}
\docType{class}
\alias{CES-class}
%\alias{calcMargins,CES-method}
%\alias{calcShares,CES-method}
%\alias{calcSlopes,CES-method}
%%\alias{CV,CES-method}
%%\alias{elast,CES-method}

\title{Class \dQuote{CES}}
\description{
  The \dQuote{CES} class contains all the information needed to calibrate a CES
  demand system and perform a merger simulation analysis under the assumption that
  firms are playing a differentiated Bertrand pricing game.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by using the constructor function \code{\link{ces}}.
%%  ~~ describe objects here ~~
}
\section{Slots}{
  Let k denote the number of products produced by all firms.
  \describe{
    \item{\code{slopes}:}{A list containing the coefficient on the numeraire
      (`alpha'),  the coefficient on price (\sQuote{gamma}), and the vector of mean
      valuations (\sQuote{meanval})}
     \item{\code{priceOutside}:}{The price of the outside good. Default
      is 1.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{Logit}}, directly.
  Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Logit}}, distance 2.
  Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 3.
}
\section{Methods}{
  For all of methods containing the \sQuote{preMerger} argument, \sQuote{preMerger} takes
  on a value of TRUE or FALSE, where TRUE invokes the method using the
  pre-merger ownership structure, while FALSE invokes the method using
  the post-merger ownership structure.
  \describe{
    \item{\code{calcShares}}{\code{signature(object, preMerger
	= TRUE,revenue=FALSE)}}{ Compute either pre-merger or post-merger equilibrium
      revenue shares under the assumptions that consumer demand is CES and firms play a differentiated product
      Bertrand Nash pricing game. \sQuote{revenue} takes
      on a value of TRUE or FALSE, where TRUE calculates revenue shares,
      while FALSE calculates quantity shares.}
    \item{\code{calcSlopes}}{\code{signature(object)}}{ Uncover CES demand
      parameters. Assumes that firms are currently at equilibrium  in a
      differentiated product Bertrand Nash pricing game.}
    \item{\code{\link{CV}}}{\code{signature(object, revenueInside)}}{
      Calculates compensating variation. If \sQuote{revenueInside} is missing, then CV
    returns compensating variation as a percent of the representative
    consumer's income. If \sQuote{revenueInside} equals the
    total expenditure on all products inside the market, then CV returns
    compensating variation in levels.}
    \item{\code{\link{elast}}}{\code{signature(object, preMerger = TRUE)}}{ Computes a k x k matrix of own and cross-price elasticities.}
	 }
}


\author{Charles Taragin \email{charles.taragin@usdoj.gov}}


\examples{
showClass("CES")           # get a detailed description of the class
showMethods(classes="CES") # show all methods defined for the class
}
\keyword{classes}
