\name{Auction2ndLogitALM-class}
\Rdversion{1.0}
\docType{class}
\alias{Auction2ndLogitALM-class}
%\alias{calcShares,LogitALM-method}
%\alias{calcSlopes,LogitALM-method}
%\alias{getNestsParms,LogitALM-method}
%%\alias{CV,LogitALM-method}
%%\alias{elast,LogitALM-method}

\title{Class \dQuote{Auction2ndLogitALM}}
\description{
  The \dQuote{Auction2ndLogitALM} class contains all the information needed to
  calibrate a Logit
  demand system with unobserved outside share and perform a merger simulation analysis under the assumption that
  firms are setting offers in a 2nd-score auction.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by using the constructor function \code{\link{auction2nd.logit.alm}}.
%%  ~~ describe objects here ~~
}

\section{Extends}{
Class \code{\linkS4class{Auction2ndLogit}}, directly.
Class \code{\linkS4class{Logit}}, distance 2.
Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Logit}}, distance 3.
Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 4.
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
    
  }

}

\author{Charles Taragin \email{charles.taragin@usdoj.gov}}

\examples{
showClass("Auction2ndLogitALM")           # get a detailed description of the class
showMethods(classes="Auction2ndLogitALM") # show all methods defined for the class
}
\keyword{classes}
