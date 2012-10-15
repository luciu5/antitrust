\name{CESNests-class}
\Rdversion{1.1}
\docType{class}
\alias{CESNests-class}
%\alias{calcShares,CESNests-method}
%\alias{calcSlopes,CESNests-method}
%\alias{CV,CESNests-method}
%\alias{elast,CESNests-method}

\title{Class "CESNests"}
\description{
  The "CESNests" class contains all the information needed to
  calibrate a nested CES
  demand system and perform a merger analysis under the assumption that
  firms are playing a differentiated products Bertrand pricing game.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by using the constructor function \code{\link{ces.nests}}.
%%  ~~ describe objects here ~~
}
\section{Slots}{
  Let k denote the number of products produced by all firms.
  \describe{
 \item{\code{nests}:}{A length k vector identifying the nest that each
    product belongs to.}
  \item{\code{parmsStart}:}{A length k vector who elements equal an
    initial guess of the nesting parameter values.}
   \item{\code{slopes}:}{A list containing the coefficient on the numeraire
      (`alpha'),  the coefficient on price (`gamma'), a vector of
      nesting parameters (`sigma') and the vector of mean
      valuations (`meanval')}
  }
}
\section{Extends}{
Class \code{\linkS4class{CES}}, directly.
Class \code{\linkS4class{Logit}}, by class "CES", distance 2.
Class \code{\linkS4class{Bertrand}}, by class "CES", distance 3.
}
\section{Methods}{
  For all of methods containing the `preMerger' argument, `preMerger' takes
  on a value of TRUE or FALSE, where TRUE invokes the method using the
  pre-merger ownership structure, while FALSE invokes the method using
  the post-merger ownership structure.
  \describe{
    \item{\code{calcShares}}{\code{signature(object, preMerger
	= TRUE, revenue = FALSE)}}{ Compute either pre-merger or post-merger equilibrium
      revenue shares under the assumptions that consumer demand is nested CES and firms play a differentiated product
      Bertrand Nash pricing game. `revenue' takes
      on a value of TRUE or FALSE, where TRUE calculates revenue shares,
      while FALSE calculates quantity shares.}
    \item{\code{calcSlopes}}{\code{signature(object)}}{ Uncover nested CES demand
      parameters. Assumes that firms are currently at equilibrium  in a
      differentiated product Bertrand Nash pricing game.}
    \item{\code{\link{CV}}}{\code{signature(object, revenueInside)}}{
      Calculates compensating variation. If `revenueInside' is missing, then CV
    returns compensating variation as a percent of the representative
    consumer's income. If `revenueInside' equals the
    total expenditure on all products inside the market, then CV returns
    compensating variation in levels.}
    \item{\code{\link{elast}}}{\code{signature(object, preMerger
	= TRUE)}}{ Computes a k x k matrix of own and cross-price
      elasticities.}
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
showClass("CESNests")
}
\keyword{classes}
