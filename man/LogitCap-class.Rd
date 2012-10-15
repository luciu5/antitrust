\name{LogitCap-class}
\Rdversion{1.1}
\docType{class}
\alias{LogitCap-class}
%\alias{calcShares,LogitCap-method}
%\alias{calcSlopes,LogitCap-method}
%\alias{getNestsParms,LogitCap-method}
%%\alias{CV,LogitCap-method}
%%\alias{elast,LogitCap-method}

\title{Class "LogitCap"}
\description{
  The "LogitCap" class contains all the information needed to
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
  Let k denote the number of products produced by all firms.
  \describe{
    \item{\code{mktSize}:}{A vector of length 1 equal to the number of
      consumers in the market. This count should include the
      number of consumers who purchase the outside option (if specified).}
    \item{\code{capacities}:}{A length k vector whose elements equal
      product capacities.}
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
    \item{\code{calcPrices}}{\code{signature(object, preMerger =
	TRUE)}}{
      Compute either pre-merger or post-merger equilibrium
      shares under the assumptions that consumer demand is Logit and firms play a differentiated product
      Bertrand Nash pricing game with capacity constraints.}
    \item{\code{calcMargins}}{\code{signature(object, preMerger =
	TRUE)}}{
      Computes equilibrium product margins assuming that firms are playing a
      Nash-Bertrand pricing game with capacity constraints. Note that
      margins for capacity constrained firms are not identified from the
      firm's first-order conditions, and so must be supplied by the user.
    }
    \item{\code{calcSlopes}}{\code{signature(object)}}{ Uncover Logit demand
      parameters. Assumes that firms are currently at equilibrium  in a
      differentiated product Bertrand Nash pricing game with capacity constraints.}
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
showClass("LogitCap")
}
\keyword{classes}
