\name{LogitNestsALM-class}
\Rdversion{1.1}
\docType{class}
\alias{LogitNestsALM-class}
%\alias{calcShares,LogitNests-method}
%\alias{calcSlopes,LogitNests-method}
%\alias{getNestsParms,LogitNests-method}
%%\alias{CV,LogitNests-method}
%%\alias{elast,LogitNests-method}

\title{Class \dQuote{LogitNestsALM}}
\description{
  The \dQuote{LogitNestsALM} class contains all the information needed to
  calibrate a nested Logit
  demand system under the assumption that the share of the outside
  product is not known. Once the model parameters have been calibrated,
  methods exist that perform a merger analysis under the assumption that
  firms are playing a differentiated products Bertrand pricing game.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by using the constructor function \code{\link{logit.nests.alm}}.
%%  ~~ describe objects here ~~
}

\section{Extends}{
  Class \code{\linkS4class{LogitNests}}, directly.
  Class \code{\linkS4class{Logit}}, by class \code{\linkS4class{LogitNests}}, distance 2.
  Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Logit}}, distance 3.
  Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 4.

}
\section{Methods}{
    For all of methods containing the \sQuote{preMerger} argument, \sQuote{preMerger} takes
  on a value of TRUE or FALSE, where TRUE invokes the method using the
  pre-merger ownership structure, while FALSE invokes the method using
  the post-merger ownership structure.
  \describe{
     \item{\code{calcSlopes}}{\code{signature(object)}}{ Uncover nested Logit demand
      parameters. Assumes that firms are currently at equilibrium  in a
      differentiated product Bertrand Nash pricing game.}
    \item{\code{\link{elast}}}{\code{signature(object, preMerger
	= TRUE)}}{ Computes a k x k matrix of own and cross-price
      elasticities.}
  }

}

\author{Charles Taragin \email{charles.taragin@usdoj.gov}}

\examples{
showClass("LogitNestsALM")           # get a detailed description of the class
showMethods(classes="LogitNestsALM") # show all methods defined for the class
}
\keyword{classes}
