\name{Stackelberg-class}
\Rdversion{1.0}
\docType{class}
\alias{Stackelberg-class}



\title{Class \dQuote{Stackelberg}}

\description{
  The \dQuote{Stackelberg} class is a building block used to create other classes in this package. As such, it is most likely to be useful for developers who wish to code their own merger calibration/simulation routines.
}

\section{Objects from the Class}{

Objects can be created by calls of the form \code{new("Stackelberg", ...)}.

}

\section{Slots}{
  Let k denote the number of products and n denote the number of plants.
  \describe{
    \item{\code{dmcfunPre}:}{A length n list whose elements equal a function that calculates the derivative of a plant's pre-merger marginal cost with respect to that plant's output. }
    \item{\code{dmcfunPost}:}{A length n list whose elements equal a function that calculates the derivative of a plant's post-merger marginal cost with respect to that plant's output. }
    \item{\code{isLeaderPre}:}{An n x k logical matrix qhose elements are TRUE if a plant produces a product pre-merger and FALSE otherwise.}
    \item{\code{isLeaderPost}:}{An n x k logical matrix qhose elements are TRUE if a plant produces a product post-merger and FALSE otherwise.}
}}


 \section{Extends}{
   Class \code{\linkS4class{Cournot}}, directly.
   Class \code{\linkS4class{Bertrand}}, by class \code{\linkS4class{Cournot}}, distance 2.
      Class \code{\linkS4class{Antitrust}}, by class \code{\linkS4class{Bertrand}}, distance 3.
 }


\author{Charles Taragin \email{charles.taragin@usdoj.gov}}


\examples{
showClass("Stackelberg")           # get a detailed description of the class
showMethods(classes="Stackelberg") # show all methods defined for the class
}
\keyword{classes}
