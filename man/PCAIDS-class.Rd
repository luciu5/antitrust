\name{PCAIDS-class}
\Rdversion{1.1}
\docType{class}
\alias{PCAIDS-class}
%\alias{calcSlopes,PCAIDS-method}


\title{Class \dQuote{PCAIDS}}
\description{
  The \dQuote{PCAIDS} class contains all the information needed to
  calibrate a PCAIDS
  demand system and perform a merger analysis under the assumption that
  firms are playing a differentiated Bertrand products pricing game.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}

\section{Objects from the Class}{
Objects can be created by using the constructor \code{\link{pcaids}}.
}

\section{Slots}{
  Let k denote the number of products produced by all firms.
  \describe{
    \item{\code{knownElast}:}{A negative number equal to the pre-merger own-price
      elasticity for any of the k products.}
    \item{\code{knownElastIndex}:}{An integer equal to the position of the
      `knownElast' product in the \sQuote{shares} vector.}

  }
}
\section{Extends}{
  Class \code{\linkS4class{AIDS}}, directly.
  Class \code{\linkS4class{Linear}}, by class \code{\linkS4class{AIDS}}, distance 2.
  Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Linear}}, distance 3.
  Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 4.
}
\section{Methods}{

  \describe{
    \item{\code{calcSlopes}}{\code{signature(object)}}{ Uncover nested CES demand
      parameters. Assumes that firms are currently at equilibrium  in a
      differentiated product Bertrand Nash pricing game.}
	 }
}

\author{Charles Taragin \email{charles.taragin@usdoj.gov}}


\examples{
showClass("PCAIDS")           # get a detailed description of the class
showMethods(classes="PCAIDS") # show all methods defined for the class
}
\keyword{classes}
