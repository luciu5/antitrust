\name{LogitNests-class}
\Rdversion{1.1}
\docType{class}
\alias{LogitNests-class}
%\alias{calcShares,LogitNests-method}
%\alias{calcSlopes,LogitNests-method}
%\alias{getNestsParms,LogitNests-method}
%%\alias{CV,LogitNests-method}
%%\alias{elast,LogitNests-method}

\title{Class \dQuote{LogitNests}}
\description{
  The \dQuote{LogitNests} class contains all the information needed to
  calibrate a nested Logit
  demand system and perform a merger simulation analysis under the assumption that
  firms are playing a differentiated products Bertrand pricing game.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by using the constructor function \code{\link{logit.nests}}.
%%  ~~ describe objects here ~~
}
\section{Slots}{
  Let k denote the number of products produced by all firms.
  \describe{
    \item{\code{nests}:}{A length k vector identifying the nest that each
    product belongs to.}
    \item{\code{parmsStart}:}{A length k vector who elements equal an
      initial guess of the nesting parameter values.}
    \item{\code{constraint}:}{A length 1 logical vector that equals TRUE
      if all nesting parameters are constrained to equal the same value
      and FALSE otherwise. Default is TRUE.}
  }
}
\section{Extends}{
Class \code{\linkS4class{Logit}}, directly.
Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Logit}}, distance 2.
}
\section{Methods}{
    For all of methods containing the \sQuote{preMerger} argument, \sQuote{preMerger} takes
  on a value of TRUE or FALSE, where TRUE invokes the method using the
  pre-merger ownership structure, while FALSE invokes the method using
  the post-merger ownership structure.
  \describe{
    \item{\code{calcShares}}{\code{signature(object, preMerger
	 = TRUE,revenue = FALSE)}}{ Compute either pre-merger or post-merger equilibrium
       shares under the assumptions that consumer demand is Logit and firms play a differentiated product
       Bertrand Nash pricing game. \sQuote{revenue} takes
       on a value of TRUE or FALSE, where TRUE calculates revenue shares,
       while FALSE calculates quantity shares.}
     \item{\code{calcSlopes}}{\code{signature(object)}}{ Uncover nested Logit demand
      parameters. Assumes that firms are currently at equilibrium  in a
      differentiated product Bertrand Nash pricing game.}
    \item{\code{\link{CV}}}{\code{signature(object)}}{ Calculate the amount of money a representative
      consumer would need to be paid to be just as well off as they were before the merger.}
    \item{\code{\link{elast}}}{\code{signature(object, preMerger
	= TRUE)}}{ Computes a k x k matrix of own and cross-price
      elasticities.}
  }

}

\author{Charles Taragin \email{charles.taragin@usdoj.gov}}

\examples{
showClass("LogitNests")           # get a detailed description of the class
showMethods(classes="LogitNests") # show all methods defined for the class
}
\keyword{classes}
