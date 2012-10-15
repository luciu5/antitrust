\name{LogitALM-class}
\Rdversion{1.1}
\docType{class}
\alias{LogitALM-class}
%\alias{calcShares,LogitALM-method}
%\alias{calcSlopes,LogitALM-method}
%\alias{getNestsParms,LogitALM-method}
%%\alias{CV,LogitALM-method}
%%\alias{elast,LogitALM-method}

\title{Class "LogitALM"}
\description{
  The "LogitALM" class contains all the information needed to
  calibrate a Logit
  demand system and perform a merger analysis under the assumption that
  firms are playing a differentiated products Bertrand pricing game with
  capacity constraints.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by using the constructor function \code{\link{logit.cap}}.
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
Class \code{\linkS4class{Bertrand}}, by class "Logit", distance 2.
}
\section{Methods}{
    For all of methods containing the `preMerger' argument, `preMerger' takes
  on a value of TRUE or FALSE, where TRUE invokes the method using the
  pre-merger ownership structure, while FALSE invokes the method using
  the post-merger ownership structure.
  \describe{
    \item{\code{calcSlopes}}{\code{signature(object)}}{ Uncover Logit
      ALM demand parameters. Assumes that firms are currently at
      equilibrium  in a differentiated product Bertrand Nash pricing game with capacity constraints.}
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
showClass("LogitALM")
}
\keyword{classes}
