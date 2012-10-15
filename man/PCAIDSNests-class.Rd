\name{PCAIDSNests-class}
\Rdversion{1.1}
\docType{class}
\alias{PCAIDSNests-class}
%\alias{calcSlopes,PCAIDSNests-method}
%\alias{getNestsParms,PCAIDSNests-method}


\title{Class "PCAIDSNests"}
\description{
  The "PCAIDSNests" class contains all the information needed to
  calibrate a nested PCAIDS
  demand system and perform a merger analysis under the assumption that
  firms are playing a differentiated Bertrand products pricing game.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}

\section{Objects from the Class}{
Objects can be created by using the constructor \code{\link{pcaids.nests}}.
}

\section{Slots}{
  Let k denote the number of products produced by all firms.
  \describe{
    \item{\code{nests}:}{A length k vector identifying which nest a product
  belongs to.}
    \item{\code{nestsParms}:}{A length k vector containing nesting parameters.}

  }
}
\section{Extends}{
  Class \code{\linkS4class{PCAIDS}}, directly.
  Class \code{\linkS4class{AIDS}}, by class "PCAIDS", distance 2.
  Class \code{\linkS4class{Linear}}, by class "AIDS", distance 3.
  Class \code{\linkS4class{Bertrand}}, by class "Linear", distance 4.
}
\section{Methods}{

  \describe{
    \item{\code{calcSlopes}}{\code{signature(object)}}{ Uncover nested CES demand
      parameters. Assumes that firms are currently at equilibrium  in a
      differentiated product Bertrand Nash pricing game.}
    \item{\code{getNestsParms}}{\code{signature(object)}}{
      Returns a matrix containing the calibrated nesting parameters. }
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
showClass("PCAIDSNests")
}
\keyword{classes}
