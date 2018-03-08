\name{Auction2ndLogit-class}
\Rdversion{1.0}
\docType{class}
\alias{Auction2ndLogit-class}
%\alias{calcShares,LogitALM-method}
%\alias{calcSlopes,LogitALM-method}
%\alias{getNestsParms,LogitALM-method}
%%\alias{CV,LogitALM-method}
%%\alias{elast,LogitALM-method}

\title{Class \dQuote{Auction2ndLogit}}
\description{
  The \dQuote{Auction2ndLogit} class contains all the information needed to
  calibrate a Logit
  demand system and perform a merger simulation analysis under the assumption that
  firms are setting offers in a 2nd-score auction.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by using the constructor function \code{\link{auction2nd.logit}}.
%%  ~~ describe objects here ~~
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
    \item{\code{calcSlopes}}{\code{signature(object)}}{ Uncover Logit
      ALM demand parameters. Assumes that firms are currently at
      equilibrium  in a 2nd-score auction.}
    \item{\code{calcMargins}}{\code{signature(object)}}{ Compute model margins.}
    \item{\code{calcMC}}{\code{signature(object)}}{ Compute constant marginal costs impied by the model.}
    \item{\code{calcShares}}{\code{signature(object)}}{ Compute logit shares using cost estimates.}
    \item{\code{calcSlopes}}{\code{signature(object)}}{ Uncover Logit
      ALM demand parameters. Assumes that firms are currently at
      equilibrium  in a 2nd-score auction.}
    \item{\code{cmcr}}{\code{signature(object)}}{ Compensating marginal cost reduction is not sensible in a 2nd-score auction. Method returns error.}
        \item{\code{upp}}{\code{signature(object)}}{ upward pricing pressure is not sensible in a 2nd-score auction. Method returns error.}
    \item{\code{CV}}{\code{signature(object)}}{ Compensating variation is simply the weighted average of the price changes.}
        \item{\code{calcPricesHypoMon}}{\code{signature(object)}}{Computes prices for a subset of firms under the control of a hypothetical monopolist playing a 2nd score auction.}
      
  }

}

\author{Charles Taragin \email{charles.taragin@usdoj.gov}}

\examples{
showClass("Auction2ndLogit")           # get a detailed description of the class
showMethods(classes="Auction2ndLogit") # show all methods defined for the class
}
\keyword{classes}
