\name{PCAIDSNests-class}
\Rdversion{1.1}
\docType{class}
\alias{PCAIDSNests-class}
%\alias{calcSlopes,PCAIDSNests-method}
\alias{getNestsParms}
\alias{getNestsParms,PCAIDSNests-method}


\title{Class \dQuote{PCAIDSNests}}
\description{
  The \dQuote{PCAIDSNests} class contains all the information needed to
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
  Class \code{\linkS4class{AIDS}}, by class \code{\linkS4class{PCAIDS}}, distance 2.
  Class \code{\linkS4class{Linear}}, by class \code{\linkS4class{AIDS}}, distance 3.
  Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Linear}}, distance 4.
  Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 5.
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

\author{Charles Taragin \email{charles.taragin@usdoj.gov}}



\examples{
showClass("PCAIDSNests")           # get a detailed description of the class
showMethods(classes="PCAIDSNests") # show all methods defined for the class
}
\keyword{classes}
