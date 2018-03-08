\name{CESALM-class}
\Rdversion{1.0}
\docType{class}
\alias{CESALM-class}
%\alias{calcShares,CESALM-method}
%\alias{calcSlopes,CESALM-method}
%\alias{getNestsParms,CESALM-method}
%%\alias{CV,CESALM-method}
%%\alias{elast,CESALM-method}

\title{Class \dQuote{CESALM}}
\description{
  The \dQuote{CESALM} class contains all the information needed to
  calibrate a CES
  demand system and perform a merger simulation analysis under the assumption that
  firms are playing a differentiated products Bertrand pricing game with
  unknown market elasticity.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by using the constructor function \code{\link{ces.alm}}.
%%  ~~ describe objects here ~~
}
\section{Slots}{
  \describe{
     \item{\code{parmsStart}:}{A length 2 vector whose first element equals an
      initial guess of the price coefficient and whose second
      element equals an initial guess of the outside share. The price
      coefficient's initial value must be negative and the outside share's
      initial value must be between 0 and 1 .}
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
    \item{\code{calcSlopes}}{\code{signature(object)}}{ Uncover CES
      ALM demand parameters. Assumes that firms are currently at
      equilibrium  in a differentiated product Bertrand Nash pricing game with capacity constraints.}
  }

}

\author{Charles Taragin \email{charles.taragin@usdoj.gov}}

\examples{
showClass("CESALM")           # get a detailed description of the class
showMethods(classes="CESALM") # show all methods defined for the class
}
\keyword{classes}
