\name{Bertrand-class}
\Rdversion{1.1}
\docType{class}
\alias{Bertrand-class}


\title{Class \dQuote{Bertrand}}
\description{
  The \dQuote{Bertrand} class is a building block used to create other classes
  in this package. As such, it is most likely to be useful for developers
  who wish to code their own merger calibration/simulation routines.


}
\section{Objects from the Class}{

Objects can be created by calls of the form \code{new("Bertrand", ...)}.

}
\section{Slots}{
  Let k denote the number of products produced by all firms.
  \describe{
    \item{\code{shares}:}{A length k vector containing observed
      output. Depending upon the model, output will be measured in
      units sold, quantity shares, or revenue shares.}
    \item{\code{mcDelta}:}{A length k vector where each element equals the
    proportional change in a product's marginal costs due to
    the merger.}
    \item{\code{slopes}:}{A k x (k+1) matrix of linear demand intercepts and slope coefficients}
    \item{subset}{A vector of length k where each element equals TRUE if
     the product indexed by that element should be included in the
     post-merger simulation and FALSE if it should be excluded.}
  }
}
 \section{Extends}{
   Class \code{\linkS4class{Antitrust}}, directly.
 }
\section{Methods}{
  Many of the methods described below contain a \sQuote{preMerger} and
  \sQuote{revenue} argument. The \sQuote{preMerger} takes
  on a value of TRUE or FALSE, where TRUE invokes the method using the
  pre-merger values, while FALSE invokes the method using
  the post-merger ownership structure. The \sQuote{revenue} argument also takes
  on a value of TRUE or FALSE, where TRUE invokes the method using
  revenues, while FALSE invokes the method using quantities
  \describe{
    \item{\code{calcMC}}{\code{signature(object,preMerger=TRUE)}}{ Calculates (constant) marginal cost for each
      product. For those classes that do not require prices, returns a
      length-k vector of NAs when prices are not supplied.}
    \item{\code{calcMargins}}{\code{signature(object, preMerger
	= TRUE)}}{ Compute either pre-merger or post-merger equilibrium
      margins under the assumption that firms play a differentiated product
      Bertrand Nash pricing game.}
    \item{\code{\link{cmcr}}}{\code{signature(object)}}{ Calculates compensated marginal cost reduction,
      the percentage decrease in the marginal costs of the merging
      parties' products needed to offset a post-merger price increase.}
    \item{\code{\link{HypoMonTest}}}{\code{signature(object,prodIndex,ssnip=.05,...)}}{
      \code{HypoMonTest} implements the Hypothetical
      Monopolist Test for a given \sQuote{ssnip}.}
    \item{\code{\link{calcPriceDeltaHypoMon}}}{\code{signature(object,prodIndex,...)}}{
      Compute the proportional difference in product prices between the
      prices of products in \sQuote{prodIndex} (i.e. prices set by the
      Hypothetical Monopolist) and prices set in the pre-merger Bertrand
      equilibrium. \sQuote{...} may be used to pass arguments to the
      optimizer.}
    \item{\code{\link{diversionHypoMon}}}{\code{signature(object,prodIndex,...)}}{
      Calculates the matrix of revenue
    diversions between all products included in the merger simulation,
    \emph{irrespective of whether or not they are also included in
    \sQuote{prodIndex}}.
     }
    \item{\code{hhi}}{\code{signature(object, preMerger= TRUE,revenue=FALSE)}}{ Compute either the pre-merger or
      post-merger Herfindahl-Hirschman Index (HHI) under the assumption
      that firms play a differentiated product Bertrand Nash pricing game.}
    \item{\code{\link{diversion}}}{\code{signature(object, preMerger
	= TRUE)}}{ Computes a k x k matrix of diversion ratios.}
    \item{\code{summary}}{\code{signature(object,revenue=TRUE,shares=TRUE,parameters=FALSE,digits=2)}}{ Summarizes
      the effect of the merger, including price and revenue
      changes. Setting \sQuote{revenue} equal to FALSE reports quantities rather
      than revenues. Setting \sQuote{shares} to FALSE reports quantities
      rather than than shares (when possible). Setting \sQuote{parameters} equal to TRUE reports
      all demand parameters. \sQuote{digits} controls the number of significant
      digits reported in output.}
     \item{\code{plot}}{\code{signature(x=object,scale=.1}}{
      Use \code{\link{ggplot}} to plot pre- and post-merger demand, marginal cost and equilibria. \sQuote{scale} controls the amount above marginal cost and below equilbrium price that is plotted.
    }
    \item{\code{upp}}{\code{signature(object)}}{
      Calculate the Upwards Pricing Pressure (upp) index.
    }
     }}

\author{Charles Taragin \email{charles.taragin@usdoj.gov}}


\examples{
showClass("Bertrand")           # get a detailed description of the class
showMethods(classes="Bertrand") # show all methods defined for the class
}
\keyword{classes}
